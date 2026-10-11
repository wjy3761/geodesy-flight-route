"""777-300ER 시나리오의 분리된 출발·도착 프로파일과 직선/원호 접속.

공식 기준에서 가져온 수치와 기하 모델의 가정은 FlightConfig/README에 구분한다.
SID/STAR, 항공기 AFM/FMC 성능, 장애물 보호구역을 재현하는 운항 모델은 아니다.
수평 직선/원호를 공항 NEU 접평면에서 설계하고 목표 타원체고면에 법선 방향으로
올린다. 실제 ECEF 이동량과 누적거리는 올린 좌표에서 다시 계산한다.
"""

from dataclasses import dataclass, field
import math
import numpy as np

from airports import default_runways
from constants import FOOT, KNOT, NAUTICAL_MILE, STANDARD_GRAVITY, WGS84_A
from coordinates import geodetic_to_ecef, ecef_to_geodetic, neu_rotation, ecef_to_neu
from navigation import destination_azimuth, signed_angle, cumulative_chord_distance
from route import constrained_step, generate_cruise, CruiseResult


@dataclass(frozen=True)
class FlightConfig:
    """가변 비행 매개변수. 각도 degree, 길이 m, 속도 IAS knot(≈CAS), Mach 무차원.

    777-300ER, 무풍, ISA를 가정한다. VREF/VR은 중량·플랩 정보 없이는 정할 수
    없으므로 190/180/150 kt 등은 예시값이며 Boeing 인증 성능값으로 주장하지 않는다.
    """
    aircraft: str = 'Boeing 777-300ER (educational scenario)'
    cruise_height_m: float = 10000.
    cruise_step_m: float = 10000.
    terminal_step_m: float = 50.
    climb_step_m: float = 500.
    takeoff_roll_m: float = 2000.
    initial_climb_deg: float = 5.
    climb_deg: float = 3.
    descent_deg: float = 3.
    final_approach_deg: float = 3.
    first_turn_height_agl_m: float = 1500.*FOOT
    minimum_turn_height_agl_m: float = 400.*FOOT
    departure_end_clearance_m: float = 35.*FOOT
    minimum_climb_gradient: float = 200.*FOOT/NAUTICAL_MILE
    # FAA 25L ILS AMDT 14C의 FAF 1900ft MSL - 시단 표고 97.8ft.
    arrival_maneuver_agl_m: float = (1900.-97.8)*FOOT
    threshold_crossing_m: float = 54.3*FOOT
    bank_limit_deg: float = 25.
    turn_rate_limit_deg_s: float = 3.
    departure_turn_ias_kt: float = 190.
    arrival_turn_ias_kt: float = 180.
    final_ias_kt: float = 150.
    climb_ias_kt: float = 290.
    low_altitude_ias_kt: float = 250.
    cruise_mach: float = 0.84
    vertical_transition_m: float = 5000.
    maneuver_north_m: float = 20000.
    maneuver_east_m: float = -20000.
    max_cruise_iterations: int = 5000
    max_connection_iterations: int = 15

    def __post_init__(self):
        for name, value in vars(self).items():
            if isinstance(value, (int, float)) and not math.isfinite(value):
                raise ValueError(f'{name} must be finite.')
        positive = ['cruise_height_m', 'cruise_step_m', 'terminal_step_m', 'climb_step_m',
                    'takeoff_roll_m', 'vertical_transition_m', 'departure_turn_ias_kt',
                    'arrival_turn_ias_kt', 'final_ias_kt', 'climb_ias_kt',
                    'low_altitude_ias_kt', 'threshold_crossing_m', 'arrival_maneuver_agl_m',
                    'minimum_turn_height_agl_m', 'departure_end_clearance_m', 'minimum_climb_gradient']
        if any(getattr(self, key) <= 0 for key in positive):
            raise ValueError('Distances, heights and speeds must be positive.')
        if not 0 < self.bank_limit_deg <= 30 or not 0 < self.turn_rate_limit_deg_s <= 3:
            raise ValueError('Turn envelope requires bank <=30 deg and rate <=3 deg/s.')
        if not 0 < self.cruise_step_m <= 100000:
            raise ValueError('Cruise step must be at most 100 km.')
        if self.first_turn_height_agl_m < self.minimum_turn_height_agl_m:
            raise ValueError('First turn is below minimum departure turn height.')
        if any(not 0 < angle < 10 for angle in [self.initial_climb_deg, self.climb_deg,
                                               self.descent_deg, self.final_approach_deg]):
            raise ValueError('Flight path angles must lie in (0,10) degrees.')
        if min(math.tan(math.radians(self.initial_climb_deg)),
               math.tan(math.radians(self.climb_deg))) < self.minimum_climb_gradient:
            raise ValueError('Configured climb is below the general 200 ft/NM benchmark.')
        if not 0 < self.cruise_mach < 1:
            raise ValueError('This ISA model supports subsonic cruise only.')
        if any(not isinstance(value, int) or isinstance(value, bool) or value <= 0
               for value in [self.max_cruise_iterations, self.max_connection_iterations]):
            raise ValueError('Iteration limits must be positive integers.')


def isa_atmosphere(height_m):
    """고도(m)의 ISA 온도 K, 압력 Pa, 음속 m/s를 반환(0~20km, 저고도는 0으로 제한).

    기하 높이를 ISA 표준 높이로 근사한다. 11km부터 등온 성층권 식을 사용한다.
    """
    if not math.isfinite(height_m) or height_m > 20000:
        raise ValueError('ISA model supports finite heights up to 20 km.')
    h = max(0., float(height_m))
    temperature = 288.15 - 0.0065*min(h, 11000.)
    pressure = 101325.*(temperature/288.15)**(STANDARD_GRAVITY/(287.05287*0.0065))
    if h > 11000:
        pressure *= math.exp(-STANDARD_GRAVITY*(h-11000)/(287.05287*temperature))
    return temperature, pressure, math.sqrt(1.4*287.05287*temperature)


def ias_to_tas(ias_kt, height_m):
    """IAS≈CAS(knots)를 ISA 압축성 피토관 식으로 TAS(m/s)로 근사 변환.

    계기·위치오차는 무시한다. 반경 계산에서 IAS를 TAS로 오인하지 않기 위한 모델.
    """
    if not math.isfinite(ias_kt) or ias_kt <= 0 or not math.isfinite(height_m):
        raise ValueError('Airspeed must be positive and parameters finite.')
    _, pressure, sound = isa_atmosphere(height_m)
    _, p0, a0 = isa_atmosphere(0.)
    impact_pressure = p0*((1.+0.2*(ias_kt*KNOT/a0)**2)**3.5 - 1.)
    mach = math.sqrt(5.*((1.+impact_pressure/pressure)**(2./7.) - 1.))
    return mach*sound


def turn_parameters(tas_mps, bank_limit_deg=25., rate_limit_deg_s=3.):
    """수평 속도(m/s), 경사각 상한(degree), 선회율 상한(degree/s)에서 반경(m)을 구함.

    β=min(β_limit, atan(Vω_limit/g)), R=V²/(g tanβ), ω=V/R.
    수평 협조 선회 식이며 작은 상승각의 출발 선회에도 기하 근사로 사용한다.
    반환값은 radius_m, bank_deg, rate_deg_s, tas_mps 딕셔너리.
    """
    if (not all(math.isfinite(v) for v in [tas_mps, bank_limit_deg, rate_limit_deg_s]) or tas_mps <= 0 or
            not 0 < bank_limit_deg < 90 or not 0 < rate_limit_deg_s):
        raise ValueError('Invalid turn parameters.')
    beta = min(math.radians(bank_limit_deg),
               math.atan(tas_mps*math.radians(rate_limit_deg_s)/STANDARD_GRAVITY))
    radius = tas_mps**2/(STANDARD_GRAVITY*math.tan(beta))
    return dict(radius_m=radius, bank_deg=math.degrees(beta),
                rate_deg_s=math.degrees(tas_mps/radius), tas_mps=tas_mps)


class AirportPlane:
    """공항 고정 NEU 접평면 좌표(n,e,m)를 지정 타원체고면에 올리는 국소 차트."""

    def __init__(self, origin):
        self.origin = tuple(origin)
        self.xyz = geodetic_to_ecef(*origin)
        self.rotation = neu_rotation(*origin[:2])

    def to_geodetic(self, xy, height):
        """평면 [East,North](m), 목표 높이(m)를 측지좌표(degree,degree,m)로 변환."""
        e, n = xy
        if math.hypot(e, n) > 500000:
            raise ValueError('Airport plane model is limited to a 500 km neighborhood.')
        up = height-self.origin[2]-(n*n+e*e)/(2*WGS84_A)
        for _ in range(15):
            xyz = self.xyz+self.rotation.T@np.array([n, e, up])
            geo = ecef_to_geodetic(*xyz)
            residual = geo[2]-height
            if abs(residual) < 1e-8:
                return geo
            # ∂h/∂u = 새 위치 법선 · 공항 원점의 Up. 고도를 임의 보정하지 않는다.
            derivative = neu_rotation(*geo[:2])[2] @ self.rotation[2]
            if derivative < 0.9:
                raise ArithmeticError('Airport height lift is outside its stable neighborhood.')
            up -= residual/derivative
        raise ArithmeticError('Airport plane height lift did not converge.')

    def to_xy(self, point):
        """측지좌표(degree,degree,m)를 공항 차트 [East,North](m)로 투영."""
        n, e, _ = ecef_to_neu(*geodetic_to_ecef(*point), *self.origin)
        return np.array([e, n])


def _theta(azimuth):
    return math.pi/2-math.radians(azimuth)


def _unit(theta):
    return np.array([math.cos(theta), math.sin(theta)])


@dataclass
class PlaneSegment:
    """평면 직선 또는 원호. 좌표 [East,North] m, 접선각 radian, 부호 +1 좌선회."""
    start: np.ndarray
    theta: float
    length_m: float
    phase: str
    radius_m: float = math.inf
    turn_sign: int = 0

    def evaluate(self, distance):
        """시작부터 거리(m)에서 평면 위치(m)와 접선각(radian)을 반환."""
        s = float(distance)
        if self.turn_sign == 0:
            return self.start+s*_unit(self.theta), self.theta
        end_theta = self.theta+self.turn_sign*s/self.radius_m
        factor = self.radius_m/self.turn_sign
        delta = factor*np.array([math.sin(end_theta)-math.sin(self.theta),
                                  -math.cos(end_theta)+math.cos(self.theta)])
        return self.start+delta, end_theta


def tangent_connection(start, start_theta, end, end_theta, radius_m, phase='arrival_turn'):
    """두 평면 위치/접선(rad)을 같은 반경(m)의 원호–직선–원호로 접선 연결.

    네 좌우 조합의 공통 전진 접선을 직접 구성하고 가장 짧은 CSC 해를 선택한다.
    Dubins 전체 6종 최단 경로 해법은 아니며 CCC 경로를 사용하지 않는다.
    반환은 PlaneSegment 리스트. 이 알고리즘은 공항 연결 구간에서만 사용한다.
    """
    if not math.isfinite(radius_m) or radius_m <= 0:
        raise ValueError('Connection radius must be positive.')
    start, end = np.asarray(start, float), np.asarray(end, float)
    candidates = []
    left = lambda theta: np.array([-math.sin(theta), math.cos(theta)])
    for s1 in [-1, 1]:
        for s2 in [-1, 1]:
            c1 = start+s1*radius_m*left(start_theta)
            c2 = end+s2*radius_m*left(end_theta)
            delta = c2-c1
            center_distance = float(np.linalg.norm(delta))
            offset = (s1-s2)*radius_m
            if center_distance < 1e-9 or abs(offset) >= center_distance:
                continue
            tangent_theta = math.atan2(delta[1], delta[0])+math.asin(offset/center_distance)
            angle1 = ((tangent_theta-start_theta)*s1) % (2*math.pi)
            angle2 = ((end_theta-tangent_theta)*s2) % (2*math.pi)
            straight = math.sqrt(max(0., center_distance**2-offset**2))
            first = PlaneSegment(start, start_theta, radius_m*angle1, phase, radius_m, s1)
            first_end, _ = first.evaluate(first.length_m)
            middle = PlaneSegment(first_end, tangent_theta, straight, 'arrival_level')
            middle_end, _ = middle.evaluate(middle.length_m)
            last = PlaneSegment(middle_end, tangent_theta, radius_m*angle2, phase, radius_m, s2)
            last_end, last_theta = last.evaluate(last.length_m)
            if np.linalg.norm(last_end-end) > 1e-5 or abs(signed_angle(math.degrees(last_theta-end_theta))) > 1e-8:
                raise ArithmeticError('CSC tangent construction failed.')
            candidates.append([first, middle, last])
    if not candidates:
        raise ArithmeticError('No forward CSC airport connection is available.')
    return min(candidates, key=lambda segments: sum(s.length_m for s in segments))


@dataclass
class ProfileResult:
    """모듈별 측지 경유점, 구간 phase(n-1), 평면 위치/거리와 계산 메타데이터."""
    points: np.ndarray
    phases: list[str]
    speeds_mps: np.ndarray
    xy: np.ndarray
    plane_distances_m: np.ndarray
    metadata: dict = field(default_factory=dict)


def _sample_segments(plane, segments, height_fn, speed_fn, config):
    points, phases, speeds, xy_points, distances = [], [], [], [], []
    offset = 0.
    for segment in segments:
        if segment.length_m < 1e-8:
            continue
        spacing = config.terminal_step_m if segment.turn_sign or segment.phase in ['final_approach', 'landing_flare', 'takeoff_roll'] else config.climb_step_m
        count = max(1, math.ceil(segment.length_m/spacing))
        samples = np.linspace(0., segment.length_m, count+1)
        for local_s in samples[0 if not points else 1:]:
            xy, theta = segment.evaluate(local_s)
            total_s = offset+local_s
            height = height_fn(total_s, segment.phase)
            geo = plane.to_geodetic(xy, height)
            if points:
                phases.append(segment.phase)
            points.append(geo)
            speeds.append(speed_fn(total_s, segment.phase, height))
            xy_points.append(xy)
            distances.append(total_s)
        offset += segment.length_m
    return ProfileResult(np.array(points), phases, np.array(speeds), np.array(xy_points), np.array(distances))


def _turn_metadata(segment, parameters):
    return dict(parameters, phase=segment.phase,
                direction='left' if segment.turn_sign > 0 else 'right',
                angle_deg=math.degrees(segment.length_m/segment.radius_m),
                horizontal_arc_m=segment.length_m,
                model_time_s=segment.length_m/parameters['tas_mps'])


def build_departure(runway, cruise_target, config=FlightConfig()):
    """활주로 시단에서 이륙·상승·원호 선회 후 순항 목표 방향으로 진입한다.

    좌표 입력은 Runway와 목표 측지좌표. 반환 ProfileResult는 WGS84 위치이다.
    1500ft 초기 상승 후 선회하고 끝 5km에서 상승기울기를 0으로 줄인다.
    선회각은 출발 활주로 접선과 실제 순항 진입 방향을 고정점 반복으로 맞춘다.
    """
    plane = AirportPlane(runway.geodetic)
    h0 = runway.geodetic[2]
    if config.takeoff_roll_m >= runway.length_m:
        raise ValueError('Takeoff roll must end inside the runway.')
    first_gain = config.first_turn_height_agl_m
    initial_gradient = math.tan(math.radians(config.initial_climb_deg))
    climb_gradient = math.tan(math.radians(config.climb_deg))
    first_distance = config.takeoff_roll_m+first_gain/initial_gradient
    clearance = (runway.length_m-config.takeoff_roll_m)*initial_gradient
    if first_distance < runway.length_m or clearance < config.departure_end_clearance_m:
        raise ValueError('Departure end clearance/straight-before-turn condition failed.')
    remaining_gain = config.cruise_height_m-h0-first_gain
    easing = config.vertical_transition_m
    high_distance = remaining_gain/climb_gradient+easing/2
    if remaining_gain <= 0 or high_distance <= easing:
        raise ValueError('Not enough climb height for level-off transition.')
    turn = turn_parameters(ias_to_tas(config.departure_turn_ias_kt,
                                    runway.orthometric_height_m+first_gain),
                           config.bank_limit_deg, config.turn_rate_limit_deg_s)
    theta0 = _theta(runway.true_bearing_deg)
    straight_start = first_distance*_unit(theta0)
    end_theta = _theta(destination_azimuth(runway.geodetic, cruise_target))
    for _ in range(40):
        delta = math.radians(signed_angle(math.degrees(end_theta-theta0)))
        arc = PlaneSegment(straight_start, theta0, abs(delta)*turn['radius_m'],
                           'departure_turn', turn['radius_m'], 1 if delta >= 0 else -1)
        if arc.length_m > high_distance-easing:
            raise ValueError('Departure turn consumes too much climb distance.')
        arc_end, _ = arc.evaluate(arc.length_m)
        end_xy = arc_end+(high_distance-arc.length_m)*_unit(end_theta)
        end_geo = plane.to_geodetic(end_xy, config.cruise_height_m)
        forward = constrained_step(end_geo, destination_azimuth(end_geo, cruise_target), 10., config.cruise_height_m)
        vector = plane.to_xy(forward)-end_xy
        new_theta = math.atan2(vector[1], vector[0])
        difference = math.radians(signed_angle(math.degrees(new_theta-end_theta)))
        if abs(difference) < 1e-8:
            break
        end_theta += difference
    else:
        raise ArithmeticError('Departure tangent iteration failed.')
    segments = [PlaneSegment(np.zeros(2), theta0, config.takeoff_roll_m, 'takeoff_roll'),
                PlaneSegment(config.takeoff_roll_m*_unit(theta0), theta0,
                             first_distance-config.takeoff_roll_m, 'initial_climb'),
                arc, PlaneSegment(arc_end, end_theta, high_distance-arc.length_m, 'climb')]

    def height_fn(s, phase):
        if s <= config.takeoff_roll_m:
            return h0
        if s <= first_distance:
            return h0+(s-config.takeoff_roll_m)*initial_gradient
        t = s-first_distance
        if t <= high_distance-easing:
            return h0+first_gain+climb_gradient*t
        q = t-(high_distance-easing)
        return h0+first_gain+climb_gradient*(high_distance-easing+q-q*q/(2*easing))

    def speed_fn(s, phase, h):
        if phase == 'departure_turn':
            return turn['tas_mps']  # 원호 중 일정 TAS 가정; IAS는 고도에 따라 달라짐.
        agl = h-h0
        ias = config.departure_turn_ias_kt if agl <= first_gain else (config.low_altitude_ias_kt if agl < 10000*FOOT else config.climb_ias_kt)
        return ias_to_tas(ias, h-runway.geoid_m)

    result = _sample_segments(plane, segments, height_fn, speed_fn, config)
    result.metadata = dict(runway=runway.to_dict(), turn=_turn_metadata(arc, turn),
                           first_turn_height_agl_m=first_gain,
                           departure_end_clearance_m=clearance,
                           initial_climb_gradient=initial_gradient,
                           climb_gradient=climb_gradient,
                           horizontal_length_m=first_distance+high_distance,
                           initial_climb_deg=config.initial_climb_deg, climb_deg=config.climb_deg,
                           level_off_transition_m=easing,
                           nadp_note='1500 ft milestone from archived ICN NADP1 guidance; thrust/flap schedule not simulated')
    return result


def _descent_height(s, length, easing, start_h, gradient):
    # 처음/끝 easing 길이에서 기울기 0↔gradient 선형 변화: 높이는 연속 이차식.
    if s <= easing:
        drop = gradient*s*s/(2*easing)
    elif s <= length-easing:
        drop = gradient*(s-easing/2)
    else:
        q = s-(length-easing)
        drop = gradient*(length-1.5*easing+q-q*q/(2*easing))
    return start_h-drop


def build_arrival(runway, source, config=FlightConfig(), inbound_theta=None):
    """순항고도 TOD에서 하강→수평 원호 연결→3° 최종 접근→플레어를 생성.

    source는 출발 쪽 측지위치, inbound_theta는 LAX 차트의 실제 진입 접선(rad).
    하강은 선회 전에 FAA 예시 FAF 고도 1900ft MSL에서 끝나고,
    접선 CSC 연결로 활주로에 정렬한다.
    시단을 TCH에서 통과한 후 이차 플레어로 활주로 내부에서 타원체고 h0에 접지한다.
    """
    plane = AirportPlane(runway.geodetic)
    h0 = runway.geodetic[2]
    maneuver_h = h0+config.arrival_maneuver_agl_m
    gradient = math.tan(math.radians(config.descent_deg))
    easing = config.vertical_transition_m
    descent_length = (config.cruise_height_m-maneuver_h)/gradient+easing
    if descent_length <= 2*easing:
        raise ValueError('Not enough height for the descent transitions.')
    maneuver = np.array([config.maneuver_east_m, config.maneuver_north_m])
    if inbound_theta is None:
        vector = maneuver-plane.to_xy(source)
        inbound_theta = math.atan2(vector[1], vector[0])
    gate_xy = maneuver-descent_length*_unit(inbound_theta)
    gate = plane.to_geodetic(gate_xy, config.cruise_height_m)
    theta_final = _theta(runway.true_bearing_deg)
    final_gradient = math.tan(math.radians(config.final_approach_deg))
    final_length = (config.arrival_maneuver_agl_m-config.threshold_crossing_m)/final_gradient
    if final_length <= 0:
        raise ValueError('Final approach start must be above threshold crossing height.')
    final_start = -final_length*_unit(theta_final)
    turn = turn_parameters(ias_to_tas(config.arrival_turn_ias_kt,
                                    runway.orthometric_height_m+config.arrival_maneuver_agl_m),
                           config.bank_limit_deg, config.turn_rate_limit_deg_s)
    connection = tangent_connection(maneuver, inbound_theta, final_start, theta_final, turn['radius_m'])
    connection_length = sum(segment.length_m for segment in connection)
    flare_length = 2*config.threshold_crossing_m/final_gradient
    if flare_length >= runway.length_m:
        raise ValueError('Flare would touch down beyond the runway.')
    segments = [PlaneSegment(gate_xy, inbound_theta, descent_length, 'descent'), *connection,
                PlaneSegment(final_start, theta_final, final_length, 'final_approach'),
                PlaneSegment(np.zeros(2), theta_final, flare_length, 'landing_flare')]
    final_offset = descent_length+connection_length
    flare_offset = final_offset+final_length

    def height_fn(s, phase):
        if s <= descent_length:
            return _descent_height(s, descent_length, easing, config.cruise_height_m, gradient)
        if s <= final_offset:
            return maneuver_h
        if s <= flare_offset:
            return maneuver_h-(s-final_offset)*final_gradient
        fraction = min(1., (s-flare_offset)/flare_length)
        return h0+config.threshold_crossing_m*(1-fraction)**2

    def speed_fn(s, phase, h):
        if phase in ['arrival_turn', 'arrival_level']:
            return turn['tas_mps']
        if phase in ['final_approach', 'landing_flare']:
            return ias_to_tas(config.final_ias_kt, h-runway.geoid_m)
        return ias_to_tas(config.low_altitude_ias_kt if h-runway.geoid_m < 10000*FOOT else config.climb_ias_kt, h-runway.geoid_m)

    result = _sample_segments(plane, segments, height_fn, speed_fn, config)
    result.metadata = dict(runway=runway.to_dict(), gate_geodetic=gate,
                           inbound_theta_rad=inbound_theta,
                           turns=[_turn_metadata(segment, turn) for segment in connection if segment.turn_sign and segment.length_m > 1e-8],
                           descent_horizontal_m=descent_length,
                           maneuver_horizontal_m=connection_length,
                           maneuver_height_agl_m=config.arrival_maneuver_agl_m,
                           final_horizontal_m=final_length,
                           final_approach_deg=config.final_approach_deg,
                           threshold_crossing_m=config.threshold_crossing_m,
                           threshold_geodetic=runway.geodetic,
                           threshold_aircraft_geodetic=plane.to_geodetic(np.zeros(2), h0+config.threshold_crossing_m),
                           touchdown_geodetic=tuple(result.points[-1]),
                           touchdown_from_threshold_m=flare_length,
                           descent_transition_m=easing,
                           profile_note='Synthetic vector connection; no published SID/STAR or terrain clearance claim')
    return result


@dataclass
class FlightResult:
    """통합 측지 경유점(n,3), 구간 phase(n-1), 누적 ECEF 거리(m), 검증 메타데이터."""
    points: np.ndarray
    phases: list[str]
    distance_m: np.ndarray
    departure: ProfileResult
    cruise: CruiseResult
    arrival: ProfileResult
    config: FlightConfig
    connection_iterations: int


def simulate_flight(config=FlightConfig()):
    """ICN33L→LAX25L의 전체 경로를 생성하고 순항/TOD 수평 접선을 반복 접속한다.

    순항은 모든 반복에서 이동 NEU 기반으로 직접 계산한다. 도착 차트에서 실제
    마지막 순항 구간의 끝 접선을 추정하여 TOD를 갱신하며 위치변화 1mm로 수렴한다.
    """
    icn, lax = default_runways()
    plane = AirportPlane(lax.geodetic)
    inbound_theta = None
    for iteration in range(1, config.max_connection_iterations+1):
        arrival = build_arrival(lax, icn.geodetic, config, inbound_theta)
        gate = tuple(arrival.points[0])
        departure = build_departure(icn, gate, config)
        cruise = generate_cruise(tuple(departure.points[-1]), gate, config.cruise_step_m,
                                 config.cruise_height_m, config.max_cruise_iterations)
        # 실제 마지막 순항 현의 끝 방향을 도착 차트에 투영하여 접속한다.
        # 현 자체의 투영은 샘플 간격에 따른 약간의 방향 오차를 남긴다.
        vector = plane.to_xy(cruise.points[-1])-plane.to_xy(cruise.points[-2])
        theta = math.atan2(vector[1], vector[0])
        next_arrival = build_arrival(lax, tuple(departure.points[-1]), config, theta)
        new_gate = np.array(next_arrival.points[0])
        shift = np.linalg.norm(geodetic_to_ecef(*new_gate)-geodetic_to_ecef(*gate))
        if shift < 1e-3:
            break
        inbound_theta = theta
    else:
        raise ArithmeticError('Cruise/descent tangent connection did not converge.')
    # 세 구간의 공유 끝점은 한 번만 저장; 각 phase는 그 구간에 들어오는 변을 뜻한다.
    points = np.vstack([departure.points, cruise.points[1:], arrival.points[1:]])
    phases = departure.phases+['cruise']*(len(cruise.points)-1)+arrival.phases
    if len(phases) != len(points)-1:
        raise ArithmeticError('Phase/waypoint lengths disagree.')
    return FlightResult(points, phases, cumulative_chord_distance(points), departure,
                        cruise, arrival, config, iteration)
