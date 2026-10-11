"""ECEF 직선의 NEU 수평 투영을 이용한 항법. 측지선 해법은 사용하지 않는다."""

import math
import numpy as np

from coordinates import geodetic_to_ecef, ecef_to_neu


def signed_angle(angle):
    """각도 차이(degree)를 [-180,180)로 정규화한다."""
    return (float(angle) + 180.0) % 360.0 - 180.0


def chord_distance(start, end):
    """두 측지좌표 (위도 degree, 경도 degree, 높이 m)의 ECEF 직선거리(m)."""
    return float(np.linalg.norm(geodetic_to_ecef(*end) - geodetic_to_ecef(*start)))


def destination_azimuth(current, destination):
    """현재 위치에서 목적지 ECEF 현의 수평 투영 방위각(degree, 북=0)을 반환.

    두 입력은 측지좌표 (위도 degree, 경도 degree, 높이 m).
    동일점·대척점 등 수평 벡터가 1μm 미만이면 방향 미정의 예외를 낸다.
    """
    target = geodetic_to_ecef(*destination)
    north, east, _ = ecef_to_neu(*target, *current)
    if math.hypot(north, east) < 1e-6:
        raise ValueError("Horizontal destination direction is undefined.")
    return math.degrees(math.atan2(east, north)) % 360.0


def cumulative_chord_distance(points):
    """측지좌표 배열 (n,3)에서 시작 0인 누적 ECEF 직선거리 배열(m)을 반환."""
    if len(points) == 0:
        return np.empty(0)
    xyz = np.array([geodetic_to_ecef(*point) for point in points])
    return np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(xyz, axis=0), axis=1))]

