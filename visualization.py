"""WGS84 ECEF, 타원체 메르카토르, 거리–타원체고 및 공항 상세 그림."""

from pathlib import Path
import math
import os
import struct
import zipfile
import numpy as np

# 제한된 Windows 환경에서도 글꼴 캐시가 프로젝트 내부에 생성되도록 한다.
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parent/'.cache'/'matplotlib'))
import matplotlib.pyplot as plt

from constants import WGS84_A, WGS84_B, WGS84_E2
from coordinates import geodetic_to_ecef

DATA = Path(__file__).resolve().parent/'data'
PHASE_STYLE = {
    'takeoff_roll': ('#8b5a2b', 'Takeoff roll'),
    'initial_climb': ('#e88d27', 'Initial climb'),
    'departure_turn': ('#cd5635', 'Departure turn'),
    'climb': ('#e8a522', 'Climb'),
    'cruise': ('#167ba8', 'Cruise'),
    'descent': ('#8c57aa', 'Descent'),
    'arrival_turn': ('#bd4480', 'Arrival turns'),
    'arrival_level': ('#b57a9d', 'Arrival level'),
    'final_approach': ('#37825a', 'Final approach'),
    'landing_flare': ('#18583b', 'Landing flare'),
}


def mercator(latitude, longitude, center_lon=0.):
    """위도·경도(degree)의 WGS84 타원체 메르카토르 x,y(m)를 반환.

    경도는 중앙경선 주위 [-180,180)로 정규화. 위도는 그림용 ±85°로 제한한다.
    y=a ln[tan(π/4+φ/2)((1-e sinφ)/(1+e sinφ))^(e/2)].
    """
    phi = np.radians(np.clip(np.asarray(latitude, float), -85., 85.))
    lam = np.radians((np.asarray(longitude, float)-center_lon+180.) % 360.-180.)
    eccentricity = math.sqrt(WGS84_E2)
    sin_phi = np.sin(phi)
    y = WGS84_A*(np.log(np.tan(math.pi/4+phi/2))
                 + eccentricity/2*np.log((1-eccentricity*sin_phi)/(1+eccentricity*sin_phi)))
    return WGS84_A*lam, y


def split_longitude_path(points, center_lon=0.):
    """(위도,경도) 배열을 지도 경계에서 보간·분할한 상대경도 배열 리스트로 반환.

    반환 열은 [위도 degree, 중앙경선 상대경도 degree]. 경계의 ±180은 그대로
    보존한다. 두 끝점만 보고 반대쪽으로 연결하지 않도록 교차점 두 개를 넣는다.
    """
    points = np.asarray(points, float)
    if len(points) == 0:
        return []
    normalized = (points[:, 1]-center_lon+180.) % 360.-180.
    result, current = [], [[points[0, 0], normalized[0]]]
    for index in range(1, len(points)):
        lat0, lon0 = points[index-1, 0], normalized[index-1]
        lat1, lon1 = points[index, 0], normalized[index]
        if abs(lon1-lon0) <= 180:
            current.append([lat1, lon1])
            continue
        unwrapped = lon1-360 if lon1-lon0 > 180 else lon1+360
        boundary = -180. if unwrapped < -180 else 180.
        fraction = (boundary-lon0)/(unwrapped-lon0)
        latitude = lat0+fraction*(lat1-lat0)
        current.append([latitude, boundary])
        result.append(np.array(current))
        current = [[latitude, -boundary], [lat1, lon1]]
    result.append(np.array(current))
    return result


def coastline_paths():
    """Natural Earth 110m SHP의 PolyLine을 표준 struct로 읽어 [위도,경도] 리스트 반환."""
    path = DATA/'ne_110m_coastline.zip'
    if not path.exists():
        raise FileNotFoundError('Coastline missing; run python tools/fetch_assets.py')
    with zipfile.ZipFile(path) as archive:
        name = next(name for name in archive.namelist() if name.endswith('.shp'))
        content = archive.read(name)
    offset, paths = 100, []
    while offset < len(content):
        _, words = struct.unpack_from('>2i', content, offset)
        record = content[offset+8:offset+8+words*2]
        offset += 8+words*2
        shape = struct.unpack_from('<i', record)[0]
        if shape == 0:
            continue
        if shape != 3:
            raise ValueError('Expected coastline PolyLine shape type 3.')
        parts, count = struct.unpack_from('<2i', record, 36)
        starts = list(struct.unpack_from(f'<{parts}i', record, 44))+[count]
        vertices = np.frombuffer(record, dtype='<f8', count=count*2, offset=44+parts*4).reshape(-1, 2)
        for first, last in zip(starts[:-1], starts[1:]):
            paths.append(vertices[first:last, ::-1].copy())
    return paths


def _runs(phases):
    first = 0
    while first < len(phases):
        last = first+1
        while last < len(phases) and phases[last] == phases[first]:
            last += 1
        yield first, last, phases[first]
        first = last


def _legend_handles(phases):
    return [plt.Line2D([0], [0], color=PHASE_STYLE[phase][0], lw=2.5,
                      label=PHASE_STYLE[phase][1]) for phase in dict.fromkeys(phases)]


def plot_3d(flight):
    """실제 ECEF 비율의 WGS84 타원체와 구간별 3차원 경로 Figure 반환."""
    fig = plt.figure(figsize=(10, 9), layout='constrained')
    ax = fig.add_subplot(projection='3d')
    lon, beta = np.meshgrid(np.linspace(-math.pi, math.pi, 96), np.linspace(-math.pi/2, math.pi/2, 48))
    x = WGS84_A*np.cos(beta)*np.cos(lon)/1000
    y = WGS84_A*np.cos(beta)*np.sin(lon)/1000
    z = WGS84_B*np.sin(beta)/1000
    ax.plot_surface(x, y, z, color='#c8dfe8', alpha=0.24, linewidth=0, shade=False)
    ax.plot_wireframe(x, y, z, rstride=6, cstride=8, color='#7594a2', alpha=0.25, linewidth=0.45)
    for coast in coastline_paths():
        xyz = np.array([geodetic_to_ecef(lat, lon, 0.) for lat, lon in coast])/1000
        ax.plot(*xyz.T, color='#567b83', alpha=0.55, lw=0.5)
    xyz = np.array([geodetic_to_ecef(*point) for point in flight.points])/1000
    for first, last, phase in _runs(flight.phases):
        ax.plot(*xyz[first:last+1].T, color=PHASE_STYLE[phase][0], lw=2.3)
    for name, point in [('ICN', xyz[0]), ('LAX', xyz[-1])]:
        ax.scatter(*point, color='#173a52', s=32)
        ax.text(*point, f'  {name}', fontsize=10)
    radius = WGS84_A/1000*1.05
    ax.set(xlim=(-radius, radius), ylim=(-radius, radius), zlim=(-radius, radius),
           xlabel='ECEF X (km)', ylabel='ECEF Y (km)', zlabel='ECEF Z (km)',
           title='ICN → LAX | WGS84 ECEF\nTrue coordinate scale · no altitude exaggeration')
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=25, azim=180)
    ax.legend(handles=_legend_handles(flight.phases), loc='upper left', fontsize=8)
    return fig


def _map_axes(ax, center, coasts):
    for coast in coasts:
        for part in split_longitude_path(coast, center):
            _, y = mercator(part[:, 0], np.zeros(len(part)))
            ax.plot(WGS84_A*np.radians(part[:, 1])/1000, y/1000, color='#7c9096', lw=0.65)
    latitudes = np.array([-60, -30, 0, 30, 60, 75])
    _, ticks = mercator(latitudes, np.zeros(len(latitudes)))
    ax.set_yticks(ticks/1000, [f'{latitude}°' for latitude in latitudes])
    longitude_ticks = np.arange(-180, 181, 60)
    ax.set_xticks(WGS84_A*np.radians(longitude_ticks)/1000,
                 [f'{int((lon+center+180)%360-180)}°' for lon in longitude_ticks])
    ax.set_xlim(-math.pi*WGS84_A/1000, math.pi*WGS84_A/1000)
    ax.set_ylim(mercator(-65, 0)[1]/1000, mercator(78, 0)[1]/1000)
    ax.set_facecolor('#f1f7fa')
    ax.grid(alpha=0.20)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('Longitude · WGS84 ellipsoidal Mercator')
    ax.set_ylabel('Latitude')


def plot_map(flight):
    """날짜변경선 분할 세계지도와 태평양 중앙경선 지도 Figure 반환."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 10), layout='constrained')
    coasts = coastline_paths()
    for ax, center, title in zip(axes, [0., 180.],
            ['World map | route split at the international date line',
             'Pacific-centered map | continuous trans-Pacific view']):
        _map_axes(ax, center, coasts)
        for first, last, phase in _runs(flight.phases):
            for part in split_longitude_path(flight.points[first:last+1, :2], center):
                _, y = mercator(part[:, 0], np.zeros(len(part)))
                ax.plot(WGS84_A*np.radians(part[:, 1])/1000, y/1000,
                        color=PHASE_STYLE[phase][0], lw=2.2)
        for name, point in [('ICN', flight.points[0]), ('LAX', flight.points[-1])]:
            x, y = mercator(*point[:2], center)
            ax.scatter(x/1000, y/1000, color='#173a52', s=30, zorder=5)
            ax.annotate(name, (x/1000, y/1000), xytext=(5, -15), textcoords='offset points')
        ax.set_title(title, loc='left', fontsize=11)
    axes[1].legend(handles=_legend_handles(flight.phases), loc='lower center', ncol=5, fontsize=8)
    fig.suptitle('ICN 33L → LAX 25L · Natural Earth coastline', fontsize=15)
    return fig


def plot_altitude(flight):
    """누적 ECEF 현 거리(km)–WGS84 타원체고(m) Figure 반환, 구간별 색상 구분."""
    fig, ax = plt.subplots(figsize=(13, 4.7), layout='constrained')
    for first, last, phase in _runs(flight.phases):
        ax.plot(flight.distance_m[first:last+1]/1000, flight.points[first:last+1, 2],
                color=PHASE_STYLE[phase][0], lw=2.8)
    ax.axhline(flight.config.cruise_height_m, color='#50748a', linestyle='--', lw=0.8,
               label=f'Cruise h = {flight.config.cruise_height_m:g} m')
    ax.set(xlabel='Cumulative ECEF chord distance (km)', ylabel='WGS84 ellipsoidal height (m)',
           title='Flight distance – ellipsoidal height | 777-300ER scenario',
           ylim=(-250, flight.config.cruise_height_m*1.08))
    ax.grid(alpha=0.2)
    ax.legend(handles=_legend_handles(flight.phases), loc='lower center', ncol=5, fontsize=8)
    return fig


def plot_terminal_details(flight):
    """선회각·원호 길이와 이착륙 고도 변화를 읽을 수 있는 공항 상세 Figure 반환."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), layout='constrained')
    for column, profile in enumerate([flight.departure, flight.arrival]):
        xy = profile.xy/1000
        ax = axes[0, column]
        for first, last, phase in _runs(profile.phases):
            ax.plot(*xy[first:last+1].T, color=PHASE_STYLE[phase][0], lw=2.3)
        runway = profile.metadata['runway']
        bearing = math.radians(runway['true_bearing_deg'])
        ax.plot([0, runway['length_m']*math.sin(bearing)/1000],
                [0, runway['length_m']*math.cos(bearing)/1000], color='#243e46', lw=4)
        ax.scatter(0, 0, color='#243e46', s=20)
        if column == 0:
            ax.set_xlim(-12, 18)
            ax.set_ylim(-5, 20)
            turns = [profile.metadata['turn']]
        else:
            ax.set_xlim(-35, 25)
            ax.set_ylim(-10, 35)
            turns = profile.metadata['turns']
        description = '\n'.join(f"{turn['direction'].title()} {turn['angle_deg']:.2f}° · R {turn['radius_m']/1000:.2f} km · arc {turn['horizontal_arc_m']/1000:.2f} km" for turn in turns)
        ax.text(0.03, 0.97, description, transform=ax.transAxes, va='top', fontsize=9,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.9))
        ax.set(title=f"{runway['airport']} RWY {runway['designation']} · local NEU plan", xlabel='East (km)', ylabel='North (km)')
        ax.set_aspect('equal', adjustable='box')
        ax.grid(alpha=0.2)
        lower = axes[1, column]
        origin_h = runway['ellipsoidal_height_m']
        for first, last, phase in _runs(profile.phases):
            lower.plot(profile.plane_distances_m[first:last+1]/1000,
                       profile.points[first:last+1, 2]-origin_h,
                       color=PHASE_STYLE[phase][0], lw=2.5)
        lower.set(xlabel='Airport profile horizontal plane distance (km)',
                  ylabel='Height above threshold ellipsoidal height (m)', title='Vertical profile')
        lower.grid(alpha=0.2)
        lower.legend(handles=_legend_handles(profile.phases), fontsize=8)
    fig.suptitle('Terminal profiles · published reference values + explicit model assumptions', fontsize=14)
    return fig


def save_figures(flight, output_dir, show=False):
    """필수 세 그림과 공항 상세도를 PNG/PDF로 저장하고 선택적으로 창을 표시."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, function in [('flight_3d', plot_3d), ('flight_map', plot_map),
                           ('flight_altitude', plot_altitude), ('terminal_details', plot_terminal_details)]:
        figure = function(flight)
        for extension in ['png', 'pdf']:
            path = output/f'{name}.{extension}'
            figure.savefig(path, dpi=180)
            paths.append(str(path))
        if not show:
            plt.close(figure)
    if show:
        plt.show()
    return paths

