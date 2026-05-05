"""Load path configuration from config.yml and expose as module-level variables."""
import os
import ast
from pathlib import Path

_cfg_path = Path(__file__).parent / "config.yml"

def _strip_inline_comment(value):
    in_single = False
    in_double = False
    for idx, char in enumerate(value):
        if char == "'" and not in_double:
            in_single = not in_single
        elif char == '"' and not in_single:
            in_double = not in_double
        elif char == "#" and not in_single and not in_double:
            return value[:idx]
    return value

def _load_simple_yaml(path):
    """Parse the simple top-level key/value config used by this project."""
    cfg = {}
    for lineno, raw_line in enumerate(path.read_text().splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if ":" not in line:
            raise ValueError(f"{path.name}:{lineno}: expected 'key: value'")

        key, value = line.split(":", 1)
        key = key.strip()
        value = _strip_inline_comment(value).strip()

        if not value or value in {"null", "Null", "NULL", "~"}:
            cfg[key] = None
        elif value[0] in {"'", '"'}:
            cfg[key] = ast.literal_eval(value)
        else:
            cfg[key] = value
    return cfg

def _load_config(path):
    try:
        import yaml
    except ModuleNotFoundError:
        return _load_simple_yaml(path)

    cfg = yaml.safe_load(path.read_text())
    return cfg or {}

_cfg = _load_config(_cfg_path)

def _is_placeholder(value):
    return str(value).strip().startswith("/path/to")

def _get(key, required=True):
    val = _cfg.get(key)
    if required and (val is None or _is_placeholder(val)):
        raise ValueError(
            f"config.yml: '{key}' has not been set. "
            f"Edit figures/config.yml before running notebooks."
        )
    if val is None or _is_placeholder(val):
        return ""
    return os.path.expanduser(str(val).strip())

def _require_dir(path, key):
    if path and not os.path.isdir(path):
        raise FileNotFoundError(f"config.yml: '{key}' directory does not exist: {path}")
    return path

def _with_sep(path):
    return path.rstrip(os.sep) + os.sep if path else ""

def _natural_earth_path(resolution, category, name):
    try:
        import cartopy.io.shapereader as shpreader
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "Cartopy is required for Natural Earth map layers. "
            "Install the project environment from environment.yml."
        )

    return shpreader.natural_earth(
        resolution=resolution,
        category=category,
        name=name,
    )

def load_world():
    import geopandas as gpd
    return gpd.read_file(naturalearth_lowres_path)

def load_coastlines():
    import geopandas as gpd
    return gpd.read_file(coastlines_path)

def load_rivers():
    import geopandas as gpd
    return gpd.read_file(rivers_path)

# -----------------------------------------------------------------------
# Required — set by user in config.yml
# -----------------------------------------------------------------------
dir_in  = _require_dir(_get("dir_in"), "dir_in")
dir_out = _get("dir_out").rstrip(os.sep) + os.sep

os.makedirs(dir_out, exist_ok=True)

# -----------------------------------------------------------------------
# Auto-derived from dir_in — everything in the Zenodo package
# -----------------------------------------------------------------------
dir_flow_managed       = _with_sep(os.path.join(dir_in, "Simulated", "SimulatedFinal-HP", "flow"))
dir_flow_natural       = _with_sep(os.path.join(dir_in, "Simulated", "SimulatedFinal-YL", "flow"))
dir_usgrid             = os.path.join(dir_in, "UScells", "contributing_grids_all")
dir_calibration_results = os.path.join(dir_in, "Simulated", "calibration results")
dir_sensitivity        = os.path.join(dir_in, "WATCH-1M-Run-abcdm")
dir_yaling             = os.path.join(dir_in, "Yaling-2018-Runoff", "q_235_basins", "235_basins")
dir_grand              = os.path.join(dir_in, "GRanD_Version_1_1")

# Reference files (BasinNames235.txt, basin.csv, Grid_Areas_ID.csv live in root)
dir_ref                = dir_in
path_grdc_stations     = os.path.join(dir_in, "GRDC_stations_list.csv")
path_grdc_csv          = os.path.join(dir_in, "grdc_91basin_m3persec_1971_1990_monthly.csv")

# -----------------------------------------------------------------------
# Optional — external sources not included in the Zenodo package
# -----------------------------------------------------------------------

# SimulatedFinal-v1 baseline flow (not in Zenodo; Figures 4 & 5 use it).
# If it is not configured, use the managed flow folder available in dir_in.
_configured_dir_flow_v1 = _require_dir(_get("dir_flow_v1", required=False), "dir_flow_v1")
dir_flow_v1 = _with_sep(_configured_dir_flow_v1) if _configured_dir_flow_v1 else dir_flow_managed

# Natural Earth map layers used by Figures 3, 5, 9 and supplement maps.
# Cartopy manages these package resources/cache locations.
coastlines_path = _get("coastlines_path", required=False) or _natural_earth_path("50m", "physical", "coastline")
rivers_path = _natural_earth_path("50m", "physical", "rivers_lake_centerlines")
naturalearth_lowres_path = _natural_earth_path("110m", "physical", "land")

# Xanthos model climate forcing (Figure 6 cell 5 only)
# Files needed: gswp3_precip, ET.RS_METEO, ETscaler, penman_monteith PET
dir_xanthos_input = _require_dir(_get("dir_xanthos_input", required=False), "dir_xanthos_input")
