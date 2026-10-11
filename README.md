# 측지학 과제 2 — ICN–LAX 비행경로 시뮬레이션

WGS84 측지좌표 ↔ ECEF ↔ 이동 원점 NEU 변환을 직접 구현하고,
현재 위치에서 목적지 방향을 다시 계산하면서 10km씩 진행한다.
순항은 고도 10,000m의 하강 시작점(TOD)까지 계산하고, 별도 공항 프로파일로 착륙한다.

## 실행

현재 프로젝트에는 Python 3.13 기반 `.venv`와 필요한 라이브러리가 준비되어 있다.
VS Code에서 이 폴더를 열고 Python 인터프리터를 `.venv/Scripts/python.exe`로 선택한다.
이 컴퓨터의 PATH 기본 `python`은 QGIS용이며 일반 Python 실행이 실패하므로,
아래처럼 프로젝트 인터프리터를 직접 사용한다.

```powershell
& .\.venv\Scripts\python.exe main.py
& .\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

검증 로그를 UTF-8 파일로도 저장하려면 `& .\.venv\Scripts\python.exe tools\run_checks.py`를 실행한다.

`main.py`는 경로 생성, 독립 비교 검증, CSV/JSON/보고서, PNG/PDF 네 그림을 만든다.
VS Code의 실행 구성도 포함되어 있다. 창에서 3D 그림을 회전하려면 다음을 실행한다.

```powershell
& .\.venv\Scripts\python.exe main.py --show
```

재설치할 때는 Python 3.10 이상으로 가상환경을 생성한 후 `requirements.txt`를 설치한다.
이 컴퓨터에서 확인된 일반 Python 경로는 다음과 같다.

```powershell
& 'C:\Users\admin\AppData\Local\Programs\Python\Python313\python.exe' -m venv .venv
& .\.venv\Scripts\python.exe -m pip install -r requirements.txt
& .\.venv\Scripts\python.exe tools\fetch_assets.py
```

다운로드한 배경지도와 지오이드 자료가 `data/`에 포함되어 있어 보통 실행에는 네트워크가 필요 없다.

## 산출물

| 파일 | 내용 |
|---|---|
| `outputs/flight_3d.png`, `.pdf` | 실제 ECEF 축 비율의 WGS84 타원체와 3차원 경로 |
| `outputs/flight_map.png`, `.pdf` | 날짜변경선에서 분할한 세계지도와 태평양 중앙 지도 |
| `outputs/flight_altitude.png`, `.pdf` | 누적 ECEF 현 거리–WGS84 타원체고 |
| `outputs/terminal_details.png`, `.pdf` | 공항별 선회각·반경·길이와 고도 프로파일 |
| `outputs/flight_path.csv` | 모든 경유점의 측지좌표와 누적거리·유입 구간 종류 |
| `outputs/validation.json` | 잔차, 측지선 비교, 입력 설정·출처·데이터 해시 |
| `outputs/RESULTS.md` | 단계별 공식·검증값·이착륙 수치·한계 설명 |
| `outputs/used_config.json` | 해당 실행에 사용한 모든 매개변수 |

거리의 계산 단위는 m이다. 읽기 쉬운 그림의 가로축에는 km를 명시한다.
각도 입출력은 degree이며 내부 삼각함수만 radian이다.

## 수학적 원리

WGS84 장반경 a=6378137m, 편평률 f=1/298.257223563, e²=f(2−f)를 사용한다.
묘유선 곡률반경 ν=a/√(1−e²sin²φ)로 ECEF를 직접 계산한다.

```text
X = (ν+h) cosφ cosλ
Y = (ν+h) cosφ sinλ
Z = (ν(1−e²)+h) sinφ
```

역변환의 위도는 φ_new=atan2(Z+e²νsinφ, √(X²+Y²))를 반복해서 구한다.
높이는 극점에서 cosφ로 나누지 않는 식을 사용한다.

```text
h = √(X²+Y²) cosφ + Z sinφ − a√(1−e²sin²φ)
```

ECEF 벡터→NEU 회전행렬은 다음과 같다. 행의 순서는 North, East, Up이다.

```text
R = [−sinφ cosλ   −sinφ sinλ    cosφ]
    [−sinλ         cosλ          0   ]
    [ cosφ cosλ    cosφ sinλ     sinφ]
```

절대 위치 변환은 NEU=R(ECEF−origin_ECEF), 역변환은 ECEF=origin_ECEF+RᵀNEU이다.
`ecef_to_neu`의 입력은 절대 ECEF 위치이므로 먼저 원점을 뺀 상대벡터를 전달하면 안 된다.
벡터 회전만 필요할 때는 `neu_rotation(...) @ delta_ecef`를 사용한다.

현재 위치를 원점으로 목적지의 NEU를 계산하고 α=atan2(E,N)으로 수평 방위각을 정한다.
U는 방향 선택에 사용하지 않는다. 동일점/대척점처럼 수평 성분이 거의 없는 경우에는
미정의 방향 예외를 낸다. 극점의 지역 방향은 선택한 원점 경도에 의존한다.

직선거리 L과 고도 h_target을 동시에 충족하기 위해 수직 이동 u 한 개만 푼다.

```text
d(u) = √(L²−u²)
ΔNEU = [d cosα, d sinα, u]
p_next = p_current + Rᵀ ΔNEU
F(u) = height(ecef_to_geodetic(p_next)) − h_target = 0
```

거리 조건은 위 식 자체로 보존된다. 국소 해 근처의 부호 변화 구간을 확장해서 찾고,
이분법으로 F(u)=0을 푼다. 각 이동 후 원점·방위각을 갱신하고 잔여거리 감소를 검사한다.
마지막 현은 10km 이하로 조정한다. 순항고도를 임의로 수정해서 거리 조건을 깨뜨리지 않는다.
고도 10km에서 10km 이동의 ΔU는 약 −7.8m이며, 이는 접평면과 곡면의 차이이다.

## Boeing 777 시나리오와 비행 기준

세부 형식은 **777-300ER**, 무풍·ISA로 선택했다. 실제 활주로의 좌표와 진방위를 사용하지만
ICN 33L 출발/LAX 25L 착륙의 선택 자체는 과제 시나리오이며 당일 운항 활주로를 의미하지 않는다.

| 항목 | 기본값 | 근거와 적용 범위 |
|---|---|---|
| 출발 일반 비교 기준 | 활주로 끝 35ft, 첫 선회 전 400ft, 200ft/NM | FAA 일반 출발 기준 참고; ICN SID의 대체 규정이 아님 |
| 첫 선회 높이 | 1500ft AGL | 모델 가정; ICN 보관 NADP1의 추력 감소 권고 높이와 구분 |
| 초기/후속 상승각 | 5° / 3° | 중량별 성능을 대신하는 모델 가정, 일반 최소 기울기보다 높음 |
| 이륙 활주거리 | 2000m | 모델 가정; VR 및 실제 필요 활주거리를 계산하지 않음 |
| 선회 경사각/선회율 상한 | 25° / 3°/s | FAA 선회·홀딩 참고로 선택한 envelope; 출도착 공통 의무값이 아님 |
| 출발 선회/도착 선회 IAS | 190kt / 180kt | 모델 가정, 기종 성능표에서 확정한 값이 아님 |
| 최종 접근 IAS | 150kt | 모델 가정; 중량·플랩·바람에 따른 실제 VREF가 아님 |
| 순항 Mach | 0.84 | Boeing의 공개 777 설계 설명 참고 |
| 하강각 | 3° | 일반적인 기하 하강 모델 가정 |
| 최종 접근각 | 3.00° | FAA LAX25L ILS AMDT14C의 공개 예시값 |
| 최종 접근 진입 고도 | 1900ft MSL | 같은 FAA 예시의 FAF 고도; 실제 FAF 좌표 전체를 재현하지는 않음 |
| 시단 통과 높이 TCH | 54.3ft | 같은 FAA 절차 스냅샷 |
| 순항 타원체고 | 10000m | 과제 조건, 실제 압력고도/비행고도 FL과 다름 |

선회각과 경로 길이는 고정 표준값을 임의로 지정하지 않는다. 출발에서는 활주로 방향과
순항 진입 방향을 연결하는 원호의 각도를 계산한다. 도착에서는 원호–직선–원호의
네 좌우 조합을 구성하고, 두 끝점과 접선을 만족하는 가장 짧은 CSC 연결을 선택한다.
이 연결은 공항 주변의 단순화된 벡터 연결이며 실제 SID/STAR가 아니다.

IAS는 계기오차를 무시해 CAS로 근사하고 ISA 피토관 식으로 TAS로 변환한다.
협조 수평 선회 모델에서 β=min(β_limit,atan(Vω_limit/g)), R=V²/(g tanβ),
원호 길이=R|Δψ|이다. 출발의 작은 상승각에도 이 수평 식을 기하 근사로 적용한다.
고정 공항 NEU 평면에서 직선·원호를 설계한 뒤 목표 높이까지 Up을 Newton 반복으로
구하므로, 실제 WGS84 높이 조건을 만족한다. 곡면에 올린 원호의 반경·길이는 평면과
조금 달라지며 실제 3차원 선분 길이 합을 별도로 보고한다.

순항 진입과 하강의 양 끝에서는 5km에 걸쳐 고도 기울기를 0과 3° 사이에서 전환한다.
최종 접근에서 TCH로 시단을 통과한 뒤 h=h_runway+TCH(1−s/L_flare)²,
L_flare=2TCH/tan(3°)인 예시 플레어로 활주로 내부에서 접지한다.
초기 상승각 전환, 최종 접근 진입의 고도 기울기 전환, 속도·추력·롤인/롤아웃의
동역학은 단순화했다. 실제 항공기 성능 및 장애물·공역·관제 적합성을 주장하지 않는다.

## 공항 고도와 데이터 출처

ICN은 2020-07-30 AIP 시단 자료의 H=6.9m, 지오이드고 21.5m를 사용해 h=28.4m이다.
LAX는 FAA 시단 표고 97.8ft를 m로 바꾼 다음 NGA EGM96 15' 격자를 직접 쌍선형
보간하여 h=H+N_g로 근사한다. LAX의 h가 음수여도 MSL 표고가 음수라는 뜻은 아니다.
지역 수직 datum 및 관측 epoch 차이를 보정하는 정밀 변환은 포함하지 않았다.
공항 자료 판본/접근 자료와 근사 여부는 CSV 외에 보고서와 JSON에도 기록한다.

- [NGA WGS84](https://earth-info.nga.mil/index.php?action=wgs84&dir=wgs84)
- [ICN AIP 2020-07-30, AD2.12](https://aim.koca.go.kr/eaipPub/Package/2020-07-30/html/eAIP/KR-AD-2.RKSI-en-GB.html)
- [ICN AIP 2021-06-15 개정 비교, NADP 및 속도 기준](https://aim.koca.go.kr/eaipPub/Package/2021-06-15-AIRAC/html/eAIP/KR-AD-2.RKSI-en-GB.html?amdt=show)
- [FAA LAX AIP](https://www.faa.gov/air_traffic/publications/atpubs/aip_html/part3_ad_2.0_california.html)
- [FAA LAX25L ILS AMDT14C, 2021 공개 자료](https://aeronav.faa.gov/acifp/1E9B75E10EF14A9890E6A9D7A07325CC-LAX/CA_KLAX_ILS%20OR%20LOC%20RWY%2025L_A14C_F.pdf)
- [FAA 일반 출발 기준](https://www.faa.gov/air_traffic/publications/aim_html/chap5_section_2.html)
- [FAA 선회·홀딩 기준 참고](https://www.faa.gov/air_traffic/publications/aim_html/chap5_section_3.html)
- [Boeing 777 Mach0.84 설명](https://www.boeing.com/commercial/777/design-highlights)
- [NGA EGM96 자료의 PROJ 배포 출처](https://github.com/OSGeo/PROJ-data/blob/master/us_nga/us_nga_README.txt)
- [Natural Earth 110m 해안선, public domain](https://www.naturalearthdata.com/downloads/110m-physical-vectors/110m-coastline/)

## 설정 변경

기본 설정을 복사하거나 `config/default.json`을 편집한 새 JSON 파일을 만든다.
부분 필드만 지정해도 나머지는 기본값이 적용된다. 길이/고도는 m, 속도는 IAS kt,
각도는 degree이다. 설정을 바꾸면 결과가 해당 비행 참고 기준과 어떻게 달라지는지도
다시 확인해야 한다.

```powershell
& .\.venv\Scripts\python.exe main.py --write-default-config config\custom.json
& .\.venv\Scripts\python.exe main.py --config config\custom.json --output outputs\custom
& .\.venv\Scripts\python.exe main.py --step 5000 --output outputs\step5000
& .\.venv\Scripts\python.exe main.py --no-plots
```

## 검증과 해석

`tests/`의 21개 테스트는 알려진 축 좌표, 전역 왕복 좌표, NEU 회전,
극점/날짜변경선, 방향 특이점, 거리·고도·방향의 동시 제약, 마지막 이동,
반복 한도, 선회 물리식, CSC 접선, 표고 기준, 단조 고도와 접속을 검사한다.
`validation.py`는 GeographicLib로 독립적인 표면 측지선 기준을 만들고,
방위각·표면 투영 길이·같은 고도로 올린 참조선까지의 근사 거리를 비교한다.

기본 결과에서 NEU 방향 모델의 기준 측지선 최대 근사 이탈은 약 427m이다.
경유점 표면 투영 길이와 최단 측지선의 차이는 약 0.068m이지만, 이 작은 거리
차이가 위치별 경로 이탈이 작다는 뜻은 아니다. 5km 간격과 비교할 때도 모델 자체의
방위각 차이는 남는다. 사용자의 승인 없이 핵심 방법을 측지선 해법으로 대체하지 않았다.

일정 높이는 경유점에 적용된다. 인접 점을 직선으로 해석하면 10km 현의 중간은
일정 높이 곡면보다 약 2m 낮다. 누적 비행거리는 이러한 ECEF 현의 합이다.
3D 그림은 실제 축 비율을 유지하므로 10km 고도와 공항 선회는 지구 전체 크기에서
작게 보인다. 공항 상세도와 고도 그래프에서 해당 구간을 확인할 수 있다.

## 모듈 구성

`constants.py` → WGS84/단위, `coordinates.py` → 좌표변환,
`navigation.py` → NEU 방향과 거리, `route.py` → 고도 제약 이동 및 순항,
`geoid.py`/`airports.py` → 표고·활주로 자료, `flight_profile.py` → 이착륙·접속,
`validation.py` → 독립 비교, `visualization.py` → 그림, `main.py` → 실행·저장.
단계별 구현·검증 기록은 `DEVELOPMENT_LOG.md`, 실행별 실제 수치는 `outputs/RESULTS.md`에 있다.

