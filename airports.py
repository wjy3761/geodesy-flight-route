"""공식 공항 자료의 재현 가능한 스냅샷과 시나리오 활주로 선택.

이 자료는 최신 운항 데이터베이스가 아니다. 날씨/관제에 따른 활주로 선택은
모델링하지 않으며, 실제 좌표의 활주로 33L→25L을 과제 시나리오로 선택한다.
"""

from dataclasses import dataclass, asdict
from constants import FOOT
from geoid import geoid_undulation

ICN_SOURCE = 'https://aim.koca.go.kr/eaipPub/Package/2020-07-30/html/eAIP/KR-AD-2.RKSI-en-GB.html'
LAX_SOURCE = 'https://www.faa.gov/air_traffic/publications/atpubs/aip_html/part3_ad_2.0_california.html'
LAX_APPROACH_SOURCE = 'https://aeronav.faa.gov/acifp/1E9B75E10EF14A9890E6A9D7A07325CC-LAX/CA_KLAX_ILS%20OR%20LOC%20RWY%2025L_A14C_F.pdf'


@dataclass(frozen=True)
class Runway:
    """활주로 시단: 좌표 degree, 높이/길이 m, 진방위 degree, 자료 출처."""
    airport: str
    designation: str
    latitude: float
    longitude: float
    orthometric_height_m: float
    geoid_m: float
    true_bearing_deg: float
    length_m: float
    source: str
    source_note: str
    geoid_note: str

    @property
    def geodetic(self):
        """시단의 (위도 degree, 경도 degree, WGS84 타원체고 m)."""
        return self.latitude, self.longitude, self.orthometric_height_m + self.geoid_m

    def to_dict(self):
        """자료 출처를 포함한 직렬화 가능한 필드와 타원체고(m)를 반환."""
        return dict(asdict(self), ellipsoidal_height_m=self.geodetic[2])


def default_runways():
    """ICN 33L 출발 및 LAX 25L 착륙의 실제 시단 자료를 반환한다.

    ICN h는 AIP의 H+GUND, LAX h는 FAA 표고+EGM96 보간으로 근사한다.
    LAX의 표고와 EGM96 사이의 수직 datum 차이는 이 과제에서 모델링하지 않는다.
    """
    icn = Runway('ICN', '33L', 37+27/60+15.21/3600, 126+27/60+39.08/3600,
                 6.9, 21.5, 324.67, 3750., ICN_SOURCE,
                 'Republic of Korea AIP snapshot 2020-07-30, AD 2.12',
                 'Published threshold geoid undulation 21.5 m')
    lat, lon = 33+56/60+14.5069/3600, -(118+22/60+57.7701/3600)
    lax = Runway('LAX', '25L', lat, lon, 97.8*FOOT, geoid_undulation(lat, lon),
                 263., 11095*FOOT, LAX_SOURCE,
                 'FAA AIP web AD 2.12 snapshot accessed 2026-10-11; true bearing rounded to degrees',
                 'NGA EGM96 15-minute bilinear approximation; local vertical datum/epoch difference omitted')
    return icn, lax

