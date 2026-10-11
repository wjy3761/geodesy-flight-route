"""이동 원점 NEU, 고정 현 길이, 일정 타원체고에 의한 반복 순항 경로."""

from dataclasses import dataclass
import math
import numpy as np

from constants import WGS84_A
from coordinates import geodetic_to_ecef, ecef_to_geodetic, neu_rotation
from navigation import chord_distance, destination_azimuth


@dataclass
class CruiseResult:
    """측지 경유점(n,3), 각 위치의 목표 잔여 현 거리(m), 구간 방위각(degree)."""
    points: np.ndarray
    remaining_m: np.ndarray
    azimuths_deg: np.ndarray


def constrained_step(current, azimuth, distance=10000.0, target_height=10000.0):
    """현재 측지좌표에서 방위각(degree), 직선거리(m), 목표 타원체고(m)로 이동.

    반환은 다음 측지좌표 (degree,degree,m). d=sqrt(L²-u²)로 거리 조건을
    해석적으로 충족하고 h(p+Rᵀ[d cosα,d sinα,u])-h_target=0을 이분법으로 푼다.
    100km 이하의 국소 전진 해를 구하며, 도달 불가능한 고도는 예외 처리한다.
    """
    if not all(math.isfinite(float(v)) for v in (azimuth, distance, target_height)):
        raise ValueError("Step parameters must be finite.")
    if not 0 < distance <= 100000:
        raise ValueError("Local chord distance must lie in (0, 100000] m.")
    if abs(target_height - current[2]) >= distance:
        raise ValueError("Height change must be smaller than chord distance.")
    xyz = geodetic_to_ecef(*current)
    rotation = neu_rotation(*current[:2]).T
    alpha = math.radians(azimuth)

    def candidate(up):
        horizontal = math.sqrt(max(0.0, (distance - up) * (distance + up)))
        delta = np.array([horizontal * math.cos(alpha), horizontal * math.sin(alpha), up])
        return ecef_to_geodetic(*(xyz + rotation @ delta))

    # 작은 이동의 곡률 처짐으로 국소 해 근방을 예측하고 부호 변화 구간 확장.
    guess = target_height - current[2] - distance**2 / (2 * (WGS84_A + current[2]))
    width = max(10.0, distance**2 / WGS84_A)
    for _ in range(30):
        lower, upper = max(-distance, guess - width), min(distance, guess + width)
        f_lower = candidate(lower)[2] - target_height
        f_upper = candidate(upper)[2] - target_height
        if f_lower <= 0 <= f_upper:
            break
        width *= 2
    else:
        raise ArithmeticError("No bracket for the forward height-constrained step.")
    for _ in range(80):
        up = (lower + upper) / 2
        point = candidate(up)
        error = point[2] - target_height
        if abs(error) <= 1e-8:
            if abs(chord_distance(current, point) - distance) > 1e-4:
                raise ArithmeticError("Chord constraint failed after inverse conversion.")
            return point
        if error < 0:
            lower = up
        else:
            upper = up
    raise ArithmeticError("Height-constrained step did not converge.")


def generate_cruise(start, target, step_m=10000.0, height_m=10000.0,
                    max_iterations=5000):
    """同높이 시작/목표 측지좌표 사이의 반복 순항을 생성하여 CruiseResult 반환.

    전체 측지선을 미리 계산하지 않는다. 매 위치에서 NEU 방위각을 갱신한다.
    마지막 현 길이는 step_m 이하로 조정한다. 잔여거리 감소를 매번 검사한다.
    """
    if not math.isfinite(step_m) or not 0 < step_m <= 100000:
        raise ValueError("Cruise step must lie in (0, 100000] m.")
    if not isinstance(max_iterations, int) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer.")
    geodetic_to_ecef(*start)
    geodetic_to_ecef(*target)
    if max(abs(start[2] - height_m), abs(target[2] - height_m)) > 1e-3:
        raise ValueError("Cruise start and target must be at the configured height.")
    points, remaining, azimuths = [tuple(start)], [chord_distance(start, target)], []
    for _ in range(max_iterations):
        current, distance = points[-1], remaining[-1]
        if distance < 1e-6:
            return CruiseResult(np.array(points), np.array(remaining), np.array(azimuths))
        azimuth = destination_azimuth(current, target)
        if distance <= step_m:
            next_point = tuple(target)
        else:
            next_point = constrained_step(current, azimuth, step_m, height_m)
        next_distance = chord_distance(next_point, target)
        if next_distance >= distance - 1e-6:
            raise ArithmeticError("Cruise did not make progress toward target.")
        points.append(next_point)
        remaining.append(next_distance)
        azimuths.append(azimuth)
    raise ArithmeticError("Cruise exceeded maximum iteration count.")

