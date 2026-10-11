"""측지학 과제 2 실행, 수치 검증, CSV/JSON/보고서 및 그림 저장."""

import argparse
import csv
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import numpy as np

from flight_profile import FlightConfig, simulate_flight
from validation import flight_validation_report

ROOT = Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'.cache'/'matplotlib'))


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f'Cannot serialize {type(value)}')


def _write_report(report, output):
    cruise, departure, arrival = report['cruise'], report['departure'], report['arrival']
    integration = report['integration']
    turn_rows = []
    for label, turn in [('ICN 출발', departure['turn']),
                         *[(f'LAX 진입 {index+1}', turn) for index, turn in enumerate(arrival['turns'])]]:
        turn_rows.append(f"| {label} | {turn['direction']} | {turn['angle_deg']:.3f} | {turn['radius_m']:.3f} | {turn['horizontal_arc_m']:.3f} | {turn['realized_3d_polyline_m']:.3f} | {turn['bank_deg']:.2f} | {turn['rate_deg_s']:.3f} |")
    text = f'''# 측지학 과제 2 — 계산 및 검증 결과

시나리오: Boeing 777-300ER, 무풍·ISA, ICN 33L → LAX 25L.
활주로·접근 자료의 공개 스냅샷 수치와 명시된 모델 가정을 사용했다.
실제 SID/STAR·기종 성능표·장애물 보호구역 전체를 재현한 운항 경로는 아니다.

## 단계 2–3: WGS84 및 좌표변환

- a=6378137 m, 1/f=298.257223563, e²=f(2-f).
- ν=a/√(1−e²sin²φ), X=(ν+h)cosφcosλ, Y=(ν+h)cosφsinλ, Z=(ν(1−e²)+h)sinφ.
- 역변환은 위도 고정점 반복, 높이는 법선 투영식으로 계산했다.
- NEU는 타원체 법선 Up을 사용하고 원점을 정확히 한 번 뺐다.
- {report['coordinates']['sample_count']}점 최대 ECEF 왕복오차: {report['coordinates']['max_ecef_round_trip_m']:.3e} m.
- 최대 높이 왕복오차: {report['coordinates']['max_height_round_trip_m']:.3e} m.
- 극점은 경도 대신 ECEF 복원 오차를 검사했다. 지구 중심은 예외 처리한다.

## 단계 4–6: 방위각 및 단일 이동

- 방위각 α=atan2(E,N), 북=0°, 시계방향. U는 방향 결정에서 제외했다.
- d=√(L²−u²), ΔN=d cosα, ΔE=d sinα로 거리 L을 보존했다.
- h(p+Rᵀ[ΔN,ΔE,u])−10000=0을 부호 구간 탐색과 이분법으로 풀었다.
- 최대 순항 이동거리 잔차: {cruise['max_full_step_error_m']:.3e} m.
- 최대 순항 높이 잔차: {cruise['max_height_error_m']:.3e} m.
- 최대 수평 방위각 잔차: {cruise['max_step_azimuth_error_deg']:.3e}°.

## 단계 7–8: 순항 생성과 독립 비교

- 순항 경유점: {cruise['waypoint_count']}개, 마지막 이동: {cruise['last_step_m']:.3f} m.
- 잔여거리 엄격 감소: {cruise['strictly_decreasing_remaining_distance']}, 최종 잔여거리: {cruise['remaining_at_target_m']:.3e} m.
- 최대 위도: {cruise['max_latitude_deg']:.4f}°, 날짜변경선 교차: {cruise['dateline_crossings']}회.
- 순항 ECEF 현 거리 합: {cruise['cumulative_chord_m']/1000:.6f} km.
- WGS84 기준 표면 측지선: {cruise['reference_surface_geodesic_m']/1000:.6f} km.
- 경유점의 표면 투영 거리 합과 기준 측지선의 차이: {cruise['surface_length_excess_m']:.3f} m.
- NEU 방위각과 측지선 방위각의 최대 차이: {cruise['max_neu_vs_geodesic_azimuth_deg']:.6f}°.
- 동일 고도로 올린 기준 측지선까지 최대 근사 이탈: {cruise['max_reference_deviation_approx_m']:.3f} m.
  참조선 1km 간격 ECEF 선분에 대한 거리이며, 횡방향 이탈의 근사값이다.
- GeographicLib는 validation.py에서만 사용하며 경로 생성에는 사용하지 않았다.

## 단계 9: 선회와 고도 프로파일

β=min(25°,atan(Vω_limit/g)), ω_limit=3°/s, R=V²/(g tanβ), 선회길이=R|Δψ|.
V는 IAS를 ISA에서 CAS≈IAS로 두고 환산한 TAS이다. 선회 중에는 일정 TAS를 가정한다.
원호는 공항 평면에서 정확하며, WGS84 높이로 올린 뒤의 실제 3D 길이는 별도로 계산했다.

| 구간 | 방향 | 선회각 ° | 반경 m | 평면 원호 m | 실제 3D 선분 합 m | 경사각 ° | 선회율 °/s |
|---|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(turn_rows)}

- 출발 활주로 진방위: {departure['runway']['true_bearing_deg']:.2f}°.
- 초기 상승각: {departure['initial_climb_deg']:.1f}°, 후속 상승각: {departure['climb_deg']:.1f}° (모델 가정).
- 출발 첫 선회 높이: {departure['first_turn_height_agl_m']:.3f} m = {departure['first_turn_height_agl_m']/0.3048:.1f}ft (모델 가정).
  ICN 보관 AIP의 NADP1 1500ft 권고는 추력 감소 높이이며, 선회 개시 의무 높이로 해석하지 않았다.
- 활주로 끝 통과 높이: {departure['departure_end_clearance_m']:.3f} m; 일반 35ft 기준보다 높다.
- 초기/후속 상승기울기: {departure['initial_climb_gradient']*100:.3f}% / {departure['climb_gradient']*100:.3f}%.
  일반 출발 비교 기준 200ft/NM ≈3.292%보다 높다. 순항 진입의 고고도 level-off는 제외한다.
- 순항 진입 마지막 5km에서 상승기울기를 0으로 줄였다.
- 하강 수평 길이: {arrival['descent_horizontal_m']/1000:.3f} km; 시작·끝 5km는 기울기 전환이다.
- 접근 선회 높이: {(arrival['maneuver_height_agl_m']+arrival['runway']['orthometric_height_m'])/0.3048:.1f}ft MSL 모델값, 시단보다 {arrival['maneuver_height_agl_m']:.3f} m 높다.
- 착륙 활주로 진방위: {arrival['runway']['true_bearing_deg']:.2f}°.
- 최종 접근각: {arrival['final_approach_deg']:.2f}°, 수평 길이: {arrival['final_horizontal_m']/1000:.3f} km.
- 시단 통과 높이 TCH: {arrival['threshold_crossing_m']:.5f} m = {arrival['threshold_crossing_m']/0.3048:.2f}ft.
- 플레어 후 접지 위치: 시단에서 {arrival['touchdown_from_threshold_m']:.3f} m, 활주로 내부.
- ICN 시단 타원체고: {departure['runway']['ellipsoidal_height_m']:.3f} m = H+N_g.
- LAX 시단 타원체고: {arrival['runway']['ellipsoidal_height_m']:.3f} m.
  FAA 표고 {arrival['runway']['orthometric_height_m']:.5f} m + EGM96 지오이드고 {arrival['runway']['geoid_m']:.5f} m.
  지역 수직 datum 차이는 보정하지 않은 과제용 근사다.

## 단계 10–11: 시각화와 통합

- WGS84 3D 그림은 ECEF 축의 실제 비율을 유지하고 고도를 과장하지 않았다.
- 메르카토르는 타원체 공식을 직접 구현했다. 세계지도에서는 ±180° 교차점을 삽입·분할한다.
- 태평양 중앙 지도도 함께 제공한다. 극지 투영은 그림용 ±85°로 제한했다.
- 고도 그래프의 세로축은 압력고도가 아닌 WGS84 타원체고이다.
- 전체 경유점 {integration['waypoint_count']}개, 총 누적 ECEF 현 거리 {integration['total_chord_distance_m']/1000:.6f} km.
- 출발–순항 접속 방향차 {integration['departure_cruise_join_deg']:.6f}°.
- 순항–하강 접속 방향차 {integration['cruise_descent_join_deg']:.6f}°, 접속 반복 {integration['connection_iterations']}회.
- 수치 통합 검증: {'PASS' if report['passed'] else 'FAIL'}.

## 해석상의 한계

경유점의 순항 타원체고는 10000m이다. 경유점 사이를 직선으로 해석하면 10km 현의
중간은 약 2m 낮아진다. 최대 반복 횟수와 잔여거리 감소를 검사하지만, 대척점 방향
특이점을 모든 전역 항로에서 해결하는 알고리즘은 아니다.

777의 중량·플랩·기상에 따른 VR/VREF, 엔진 추력, 가속/감속·롤인/롤아웃은 모델링하지
않았다. 출발 190kt, 도착 선회 180kt, 최종 접근 150kt와 상승각·이륙 활주거리·플레어는
기종 성능표를 대신하는 조정 가능한 시나리오 가정이다. SID/STAR와 지형·공역 검증을
포함하지 않으므로 수치 검증 통과를 실제 운항 표준의 전체 준수로 해석해서는 안 된다.

## 자료 출처

- [WGS84 NGA](https://earth-info.nga.mil/index.php?action=wgs84&dir=wgs84)
- [ICN AIP 2020-07-30 AD 2.12]({departure['runway']['source']})
- [ICN AIP 2021-06-15 개정 비교, NADP와 속도 기준](https://aim.koca.go.kr/eaipPub/Package/2021-06-15-AIRAC/html/eAIP/KR-AD-2.RKSI-en-GB.html?amdt=show)
- [FAA LAX AIP]({arrival['runway']['source']})
- [FAA LAX 25L ILS AMDT 14C, 2021 예시: 3°, TCH 54.3ft, FAF 1900ft](https://aeronav.faa.gov/acifp/1E9B75E10EF14A9890E6A9D7A07325CC-LAX/CA_KLAX_ILS%20OR%20LOC%20RWY%2025L_A14C_F.pdf)
- [FAA 일반 출발 기준](https://www.faa.gov/air_traffic/publications/aim_html/chap5_section_2.html)
- [FAA 선회·홀딩 기준: 선회 envelope 선택의 참고, 출도착 공통 규정으로 적용하지 않음](https://www.faa.gov/air_traffic/publications/aim_html/chap5_section_3.html)
- [Boeing 777 설계: 통상 Mach 0.84](https://www.boeing.com/commercial/777/design-highlights)
- [NGA EGM96 격자 출처](https://github.com/OSGeo/PROJ-data/blob/master/us_nga/us_nga_README.txt)
- [Natural Earth 110m 해안선](https://www.naturalearthdata.com/downloads/110m-physical-vectors/110m-coastline/)
'''
    if 'half_step_sensitivity' in cruise:
        fine = cruise['half_step_sensitivity']
        text += f"\n## 이동 간격 민감도\n\n{fine['step_m']:.0f}m 이동 시 {fine['waypoint_count']}개 경유점, 누적 현 거리 변화 {fine['length_difference_m']:.3f}m.\n기존 경유점에서 세분 경로 선분까지 최대 거리 {fine['max_coarse_to_fine_polyline_m']:.3f}m.\n이 수치는 세분 선분의 곡률 처짐도 포함하므로 순수한 횡방향 차이는 아니다.\n"
    (output/'RESULTS.md').write_text(text, encoding='utf-8')


def main():
    """CLI 설정을 읽어 경로와 검증/그림 산출물을 생성한다. 오류는 종료코드에 반영."""
    parser = argparse.ArgumentParser(description='WGS84 ICN–LAX flight path simulation')
    parser.add_argument('--config', type=Path, help='FlightConfig 필드의 JSON 파일')
    parser.add_argument('--step', type=float, help='순항 ECEF 이동 간격 m')
    parser.add_argument('--output', type=Path, default=ROOT/'outputs')
    parser.add_argument('--no-plots', action='store_true')
    parser.add_argument('--show', action='store_true')
    parser.add_argument('--no-sensitivity', action='store_true')
    parser.add_argument('--write-default-config', type=Path)
    args = parser.parse_args()
    if args.write_default_config:
        args.write_default_config.write_text(json.dumps(asdict(FlightConfig()), indent=2), encoding='utf-8')
        print(f'Wrote config: {args.write_default_config}')
        return
    values = json.loads(args.config.read_text(encoding='utf-8-sig')) if args.config else {}
    if args.step is not None:
        values['cruise_step_m'] = args.step
    config = FlightConfig(**values)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    flight = simulate_flight(config)
    report = flight_validation_report(flight, sensitivity=not args.no_sensitivity)
    with (output/'flight_path.csv').open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index', 'latitude_deg', 'longitude_deg', 'ellipsoidal_height_m',
                         'cumulative_chord_m', 'incoming_phase'])
        for index, (point, distance) in enumerate(zip(flight.points, flight.distance_m)):
            writer.writerow([index, *point, distance, 'start' if index == 0 else flight.phases[index-1]])
    report['data_sha256'] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                             for path in (ROOT/'data').glob('*') if path.suffix in ['.zip', '.tif']}
    (output/'validation.json').write_text(json.dumps(report, indent=2, default=_json_default), encoding='utf-8')
    _write_report(report, output)
    (output/'used_config.json').write_text(json.dumps(asdict(config), indent=2), encoding='utf-8')
    if not args.no_plots:
        if not args.show:
            import matplotlib
            matplotlib.use('Agg')
        from visualization import save_figures
        save_figures(flight, output, args.show)
    print(f"PASS | {len(flight.points)} waypoints | {flight.distance_m[-1]/1000:.3f} km")
    print(f"Cruise max height error: {report['cruise']['max_height_error_m']:.3e} m")
    print(f"Cruise max chord error: {report['cruise']['max_full_step_error_m']:.3e} m")
    print(f"Geodesic maximum approximate deviation: {report['cruise']['max_reference_deviation_approx_m']:.3f} m")
    print(f'Results: {output}')


if __name__ == '__main__':
    main()
