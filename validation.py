"""독립 수치 검증. GeographicLib는 이 모듈에서 비교 기준으로만 사용한다."""

from dataclasses import asdict
import math
import numpy as np
from geographiclib.geodesic import Geodesic

from coordinates import geodetic_to_ecef, ecef_to_geodetic, ecef_to_neu
from navigation import chord_distance, destination_azimuth, signed_angle, cumulative_chord_distance
from route import generate_cruise


def coordinate_accuracy_report():
    """재현 가능한 전역 1000점과 극점·날짜변경선의 왕복 오차 요약(m)을 반환."""
    rng = np.random.default_rng(20261011)
    points = list(zip(rng.uniform(-90, 90, 1000), rng.uniform(-180, 180, 1000),
                      rng.uniform(-500, 1e6, 1000)))
    points += [(90, 45, 10000), (-90, -180, 0), (89.999999, 179.999999, 10000)]
    xyz_errors, height_errors = [], []
    for point in points:
        xyz = geodetic_to_ecef(*point)
        restored = ecef_to_geodetic(*xyz)
        xyz_errors.append(float(np.linalg.norm(geodetic_to_ecef(*restored)-xyz)))
        height_errors.append(abs(restored[2]-point[2]))
    return dict(sample_count=len(points), max_ecef_round_trip_m=max(xyz_errors),
                max_height_round_trip_m=max(height_errors))


def _point_to_polyline_distances(points, reference):
    # 3차원 선분에 정사영한 거리: 참조 곡선 샘플링의 현 오차도 포함한다.
    a, delta = reference[:-1], np.diff(reference, axis=0)
    denominator = np.sum(delta*delta, axis=1)
    distances = []
    for point in points:
        fraction = np.clip(np.sum((point-a)*delta, axis=1)/denominator, 0., 1.)
        projections = a+fraction[:, None]*delta
        distances.append(float(np.min(np.linalg.norm(projections-point, axis=1))))
    return np.array(distances)


def cruise_accuracy_report(cruise, step_m=10000., height_m=10000., sensitivity=True):
    """순항 거리/고도/방위각/수렴과 WGS84 측지선 비교 요약을 반환.

    비교 측지선은 동일 시작/목표의 타원체 표면 측지선이다. 이탈량은 그 측지선을
    같은 고도에 올린 1km 간격 ECEF 선분까지의 근사 거리이며 횡방향 오차의 근사값.
    surface_length_excess는 고도 효과를 제외한 경유점 지표 투영 길이 차이다.
    """
    points = cruise.points
    if len(points) < 2:
        raise ValueError('Accuracy comparison requires at least two cruise points.')
    step_lengths = np.array([chord_distance(a, b) for a, b in zip(points[:-1], points[1:])])
    direction_errors, azimuth_differences = [], []
    surface_length = 0.
    target = tuple(points[-1])
    for index, (a, b) in enumerate(zip(points[:-1], points[1:])):
        n, e, _ = ecef_to_neu(*geodetic_to_ecef(*b), *a)
        direction_errors.append(abs(signed_angle(math.degrees(math.atan2(e, n))-cruise.azimuths_deg[index])))
        reference_direction = Geodesic.WGS84.Inverse(a[0], a[1], target[0], target[1])['azi1']
        azimuth_differences.append(abs(signed_angle(cruise.azimuths_deg[index]-reference_direction)))
        surface_length += Geodesic.WGS84.Inverse(a[0], a[1], b[0], b[1])['s12']
    inverse = Geodesic.WGS84.Inverse(*points[0, :2], *points[-1, :2])
    line = Geodesic.WGS84.Line(*points[0, :2], inverse['azi1'])
    reference_geos = []
    for distance in np.linspace(0., inverse['s12'], math.ceil(inverse['s12']/1000)+1):
        point = line.Position(float(distance))
        reference_geos.append((point['lat2'], point['lon2'], height_m))
    reference_xyz = np.array([geodetic_to_ecef(*point) for point in reference_geos])
    route_xyz = np.array([geodetic_to_ecef(*point) for point in points])
    deviation = _point_to_polyline_distances(route_xyz, reference_xyz)
    report = dict(waypoint_count=len(points), full_step_count=max(0, len(step_lengths)-1),
                  max_height_error_m=float(np.max(abs(points[:, 2]-height_m))),
                  max_full_step_error_m=float(np.max(abs(step_lengths[:-1]-step_m))) if len(step_lengths) > 1 else 0.,
                  last_step_m=float(step_lengths[-1]),
                  max_step_azimuth_error_deg=max(direction_errors),
                  strictly_decreasing_remaining_distance=bool(np.all(np.diff(cruise.remaining_m)<0)),
                  remaining_at_target_m=float(cruise.remaining_m[-1]),
                  max_neu_vs_geodesic_azimuth_deg=max(azimuth_differences),
                  cumulative_chord_m=float(sum(step_lengths)),
                  surface_projected_length_m=surface_length,
                  reference_surface_geodesic_m=inverse['s12'],
                  surface_length_excess_m=surface_length-inverse['s12'],
                  max_reference_deviation_approx_m=float(max(deviation)),
                  reference_sampling_m=1000.,
                  max_latitude_deg=float(max(points[:, 0])),
                  dateline_crossings=int(np.sum(abs(np.diff(points[:, 1]))>180.)))
    if sensitivity:
        fine = generate_cruise(tuple(points[0]), target, step_m/2, height_m)
        fine_xyz = np.array([geodetic_to_ecef(*point) for point in fine.points])
        report['half_step_sensitivity'] = dict(step_m=step_m/2, waypoint_count=len(fine.points),
                length_difference_m=float(cumulative_chord_distance(fine.points)[-1]-sum(step_lengths)),
                max_coarse_to_fine_polyline_m=float(max(_point_to_polyline_distances(route_xyz, fine_xyz))),
                note='Polyline comparison includes the half-step chord sag, not purely lateral displacement.')
    return report


def _join_angle(first, second):
    origin = tuple(first[-1])
    backward = destination_azimuth(origin, tuple(first[-2]))+180.
    forward = destination_azimuth(origin, tuple(second[1]))
    return abs(signed_angle(backward-forward))


def _turn_realized_lengths(profile):
    lengths = np.diff(cumulative_chord_distance(profile.points))
    runs, active = [], []
    for index, phase in enumerate(profile.phases):
        if phase.endswith('_turn'):
            active.append(float(lengths[index]))
        elif active:
            runs.append(sum(active))
            active = []
    if active:
        runs.append(sum(active))
    return runs


def flight_validation_report(flight, sensitivity=True):
    """최종 통합 검증 결과 및 공식 수치/가정을 포함하는 직렬화 가능한 사전 반환."""
    f = flight
    cruise = cruise_accuracy_report(f.cruise, f.config.cruise_step_m,
                                    f.config.cruise_height_m, sensitivity)
    departure = dict(f.departure.metadata)
    departure_lengths = _turn_realized_lengths(f.departure)
    departure['turn'] = dict(departure['turn'], realized_3d_polyline_m=departure_lengths[0] if departure_lengths else 0.)
    arrival = dict(f.arrival.metadata)
    arrival['turns'] = [dict(turn, realized_3d_polyline_m=length)
                        for turn, length in zip(arrival['turns'], _turn_realized_lengths(f.arrival))]
    report = dict(config=asdict(f.config), coordinates=coordinate_accuracy_report(), cruise=cruise,
                  departure=departure, arrival=arrival,
                  integration=dict(waypoint_count=len(f.points), total_chord_distance_m=float(f.distance_m[-1]),
                                   departure_cruise_join_deg=_join_angle(f.departure.points, f.cruise.points),
                                   cruise_descent_join_deg=_join_angle(f.cruise.points, f.arrival.points),
                                   connection_iterations=f.connection_iterations,
                                   departure_monotonic_height=bool(np.all(np.diff(f.departure.points[:, 2])>=-1e-6)),
                                   arrival_monotonic_height=bool(np.all(np.diff(f.arrival.points[:, 2])<=1e-6))),
                  limitations=[
                      'Core cruise is ECEF-chord horizontal pursuit, not an ellipsoid geodesic solver.',
                      'Constant cruise height applies to vertices; straight chords sag between vertices.',
                      '777-300ER speed/gradient/roll/flare defaults are scenario assumptions, not AFM performance.',
                      'Published procedure numbers are snapshots; synthetic connections are not SID/STAR reproductions.',
                      'No obstacle/airspace/weather/ATC verification; no claim of operational flight compliance.',
                      'LAX orthometric-to-ellipsoidal conversion uses approximate EGM96, not a full local datum transform.',
                      'Airport arcs are exact in the local plane; lifting to WGS84 slightly changes physical radius/length.',
                      'IAS is approximated as CAS in ISA; geometric ellipsoidal height is not a pressure flight level.',
                      'Roll-in/out, engine thrust, flap dynamics and mass-specific performance are omitted.',
                      'Final approach follows a surface-height model with a simplified flare and flat runway height.'])
    checks = [cruise['max_height_error_m'] < 1e-3,
              cruise['max_full_step_error_m'] < 1e-3,
              cruise['max_step_azimuth_error_deg'] < 1e-6,
              cruise['strictly_decreasing_remaining_distance'],
              cruise['remaining_at_target_m'] < 1e-3,
              report['integration']['departure_cruise_join_deg'] < 0.1,
              report['integration']['cruise_descent_join_deg'] < 0.1,
              report['integration']['departure_monotonic_height'],
              report['integration']['arrival_monotonic_height']]
    report['passed'] = bool(all(checks))
    if not report['passed']:
        raise ArithmeticError(f'Integrated validation failed: {checks}')
    return report

