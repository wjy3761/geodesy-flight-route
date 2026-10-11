"""직접 구현한 WGS84 측지좌표, ECEF, 이동 원점 NEU 변환.

공개 함수는 단일 위치를 처리한다. 위도/경도는 degree, 높이/좌표는 m.
NEU의 Up은 지구중심 방향이 아니라 타원체의 법선 방향이다.
"""

import math
import numpy as np

from constants import WGS84_A, WGS84_B, WGS84_E2


def _finite(*values):
    if not all(math.isfinite(float(v)) for v in values):
        raise ValueError("Coordinates must be finite.")


def normalize_longitude(longitude):
    """경도(degree)를 [-180, 180)로 정규화하여 반환한다."""
    _finite(longitude)
    return (float(longitude) + 180.0) % 360.0 - 180.0


def geodetic_to_ecef(latitude, longitude, height):
    """측지위도/경도(degree), 타원체고(m)를 ECEF 배열 [X,Y,Z](m)로 변환."""
    _finite(latitude, longitude, height)
    if not -90.0 <= latitude <= 90.0:
        raise ValueError("Latitude must lie in [-90, 90] degrees.")
    phi = math.radians(latitude)
    lam = math.radians(normalize_longitude(longitude))
    sin_phi, cos_phi = math.sin(phi), math.cos(phi)
    # 묘유선 곡률반경 ν=a/sqrt(1-e²sin²φ).
    nu = WGS84_A / math.sqrt(1.0 - WGS84_E2 * sin_phi**2)
    return np.array([(nu + height) * cos_phi * math.cos(lam),
                     (nu + height) * cos_phi * math.sin(lam),
                     (nu * (1.0 - WGS84_E2) + height) * sin_phi])


def ecef_to_geodetic(x, y, z):
    """ECEF(m)를 (위도 degree, 경도 degree, 타원체고 m)로 역변환.

    위도를 고정점 반복으로 구하고, 법선에 투영하여 극점에서도 안정적으로
    높이를 계산한다. 지구 중심에서는 위치가 미정의이므로 ValueError.
    항공·지표 위치용이며 지구 깊은 내부의 비유일한 해는 지원하지 않는다.
    """
    _finite(x, y, z)
    p = math.hypot(x, y)
    if math.hypot(p, z) < 1.0:
        raise ValueError("Geodetic coordinates are undefined near Earth center.")
    if p < 1e-9:
        return math.copysign(90.0, z), 0.0, abs(z) - WGS84_B
    lam = math.atan2(y, x)
    phi = math.atan2(z, p * (1.0 - WGS84_E2))
    for _ in range(30):
        nu = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(phi)**2)
        updated = math.atan2(z + WGS84_E2 * nu * math.sin(phi), p)
        if abs(updated - phi) < 2e-15:
            phi = updated
            break
        phi = updated
    else:
        raise ArithmeticError("ECEF inverse latitude did not converge.")
    # p cosφ + z sinφ = h + a sqrt(1-e²sin²φ): cosφで除さない。
    height = (p * math.cos(phi) + z * math.sin(phi)
              - WGS84_A * math.sqrt(1.0 - WGS84_E2 * math.sin(phi)**2))
    return math.degrees(phi), normalize_longitude(math.degrees(lam)), height


def neu_rotation(origin_lat, origin_lon):
    """원점의 위도/경도(degree)로 ECEF 벡터→NEU 회전행렬(무차원)을 반환."""
    _finite(origin_lat, origin_lon)
    if not -90 <= origin_lat <= 90:
        raise ValueError("Origin latitude is out of range.")
    phi, lam = math.radians(origin_lat), math.radians(origin_lon)
    sp, cp, sl, cl = math.sin(phi), math.cos(phi), math.sin(lam), math.cos(lam)
    return np.array([[-sp * cl, -sp * sl, cp],
                     [-sl, cl, 0.0],
                     [cp * cl, cp * sl, sp]])


def ecef_to_neu(x, y, z, origin_lat, origin_lon, origin_h):
    """절대 ECEF 위치(m)를 원점 측지좌표(degree,m) 기준 NEU 위치(m)로 변환.

    입력은 상대 벡터가 아닌 절대 위치이다. 원점 ECEF를 한 번만 뺀다.
    """
    _finite(x, y, z)
    origin = geodetic_to_ecef(origin_lat, origin_lon, origin_h)
    return neu_rotation(origin_lat, origin_lon) @ (np.array([x, y, z]) - origin)


def neu_to_ecef(n, e, u, origin_lat, origin_lon, origin_h):
    """NEU 위치(m)를 원점 측지좌표(degree,m) 기준 절대 ECEF 위치(m)로 변환."""
    _finite(n, e, u)
    origin = geodetic_to_ecef(origin_lat, origin_lon, origin_h)
    return origin + neu_rotation(origin_lat, origin_lon).T @ np.array([n, e, u])

