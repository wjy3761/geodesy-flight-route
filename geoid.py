"""NGA EGM96 15' 격자의 직접 쌍선형 보간: h=H+N_g.

PROJ가 배포하는 NGA 자료만 읽으며 PyProj 좌표변환은 사용하지 않는다.
EGM96은 지역 수직 datum/관측 epoch 차이를 정밀하게 보정하는 모델은 아니다.
"""

from functools import lru_cache
from pathlib import Path
import math
import numpy as np
from PIL import Image

GRID_PATH = Path(__file__).resolve().parent / 'data' / 'us_nga_egm96_15.tif'


@lru_cache(maxsize=1)
def _load_grid():
    if not GRID_PATH.exists():
        raise FileNotFoundError("EGM96 grid missing; run python tools/fetch_assets.py")
    with Image.open(GRID_PATH) as image:
        scale = image.tag_v2[33550]
        tie = image.tag_v2[33922]
        if image.size != (1440, 721) or scale[:2] != (0.25, 0.25) or tie[3:5] != (-180., 90.):
            raise ValueError("Unexpected EGM96 grid geometry.")
        return np.asarray(image, dtype=float).copy()


def geoid_undulation(latitude, longitude):
    """위도/경도(degree)에서 EGM96 지오이드고(m)를 주기 경도 쌍선형 보간."""
    if not math.isfinite(latitude) or not math.isfinite(longitude) or not -90 <= latitude <= 90:
        raise ValueError("Invalid geoid location.")
    grid = _load_grid()
    x = ((longitude + 180) % 360) / 0.25
    y = (90 - latitude) / 0.25
    ix, iy = math.floor(x), min(math.floor(y), 719)
    fx, fy = x - ix, y - iy
    jx = (ix + 1) % 1440
    return float((1-fy) * ((1-fx)*grid[iy, ix] + fx*grid[iy, jx])
                 + fy * ((1-fx)*grid[iy+1, ix] + fx*grid[iy+1, jx]))

