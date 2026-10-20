"""Plot verified, unfiltered ASO extreme locations with the actual watershed.

The external input table contains the original paired raster samples. Nothing
is reclassified, filtered, resampled, or recomputed from model output here.
The boundary is the selected HUC10 polygon, never an outline of a raster mask.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import pyproj
from pyproj import CRS, Transformer
from pyproj.exceptions import ProjError
import shapely
from shapely.geometry import Point


HUC10 = "1402000102"
REFERENCE_FLIGHTS = (
    "2018-03-31", "2018-05-24", "2019-04-07", "2019-06-10",
    "2022-04-21", "2022-05-18", "2023-04-01", "2023-05-23",
    "2024-04-03", "2024-05-20", "2025-04-07", "2025-05-20",
    "2026-03-18", "2026-04-05",
)
THRESHOLD_COUNTS = {
    "2019-04-07": 39, "2019-06-10": 25, "2022-05-18": 3,
    "2023-04-01": 2, "2023-05-23": 8, "2024-05-20": 1,
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def software_versions():
    return {"python": sys.version, "numpy": np.__version__, "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__, "geopandas": gpd.__version__,
            "pyproj": pyproj.__version__, "shapely": shapely.__version__}


def load_locations(path):
    table = pd.read_csv(path, dtype={"date": str, "crs": str})
    required = {"date", "row", "column", "x_m", "y_m", "crs",
                "longitude", "latitude", "swe_m", "snow_depth_m"}
    if not required.issubset(table.columns):
        raise ValueError("Location table lacks required original-sample columns")
    if table.duplicated(["date", "row", "column"]).any():
        raise ValueError("Duplicate flight/raster-center records")
    fields = ["x_m", "y_m", "longitude", "latitude", "swe_m", "snow_depth_m"]
    if not np.isfinite(table[fields].to_numpy(float)).all():
        raise ValueError("Location/value table contains missing or nonfinite data")
    if not ((table.swe_m > 5) | (table.snow_depth_m > 15)).all():
        raise ValueError("Every supplied row must meet the original extreme criterion")
    if table.groupby("date").size().to_dict() != THRESHOLD_COUNTS:
        raise ValueError("The six-flight, 78-record reference diagnostic has changed")
    if table.crs.nunique() != 1:
        raise ValueError("All map coordinates must use one declared CRS")
    crs = CRS.from_user_input(table.crs.iloc[0])
    if crs != CRS.from_epsg(32613):
        raise ValueError("The reference extreme coordinates must be UTM13N EPSG:32613")
    lon, lat = Transformer.from_crs(crs, 4326, always_xy=True).transform(
        table.x_m.to_numpy(), table.y_m.to_numpy())
    if not (np.allclose(lon, table.longitude, rtol=0, atol=1e-7)
            and np.allclose(lat, table.latitude, rtol=0, atol=1e-7)):
        raise ValueError("Projected and geographic source coordinates disagree")
    return table


def recurrence_mask(table):
    """Original cross-flight proximity rule; 125 m includes its endpoint."""
    xy = table[["x_m", "y_m"]].to_numpy(float)
    different = table.date.to_numpy()[:, None] != table.date.to_numpy()[None, :]
    distance = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(axis=2))
    return ((distance <= 125) & different).any(axis=1)


def load_flights(path):
    table = pd.read_csv(path, dtype={"date": str})
    if "date" not in table or tuple(table.date) != REFERENCE_FLIGHTS:
        raise ValueError("Expected the exact 14 reference ASO flight dates in chronological order")
    return table.date.tolist()


def load_boundary(path, target_crs="EPSG:32613", huc_field="huc10"):
    source = gpd.read_file(path)
    if source.crs is None:
        raise ValueError("The watershed file must declare its source CRS")
    fields = {str(c).lower(): c for c in source.columns}
    if huc_field.lower() not in fields:
        raise ValueError("The watershed file lacks the requested HUC10 field")
    field = fields[huc_field.lower()]
    selected = source[source[field].astype(str).eq(HUC10)].copy()
    if len(selected) != 1:
        raise ValueError("Expected exactly one East River HUC10=1402000102 feature")
    geometry = selected.geometry.iloc[0]
    if geometry.geom_type not in {"Polygon", "MultiPolygon"}:
        raise ValueError("The watershed feature must be polygonal")
    if geometry.is_empty or not geometry.is_valid:
        raise ValueError("The watershed polygon is empty or invalid")
    transformed = selected.to_crs(target_crs)
    transform = Transformer.from_crs(source.crs, target_crs, always_xy=True)
    first_polygon = geometry.geoms[0] if geometry.geom_type == "MultiPolygon" else geometry
    transform.transform(*first_polygon.exterior.coords[0])
    try:
        operation = transform.get_last_used_operation()
    except ProjError:  # An already projected equivalent file uses a no-op.
        operation = transform
    metadata = {
        "huc10": HUC10,
        "attributes": {str(k): (None if pd.isna(v)
                               else v.isoformat() if isinstance(v, pd.Timestamp)
                               else v.item() if hasattr(v, "item") else v)
                       for k, v in selected.iloc[0].items() if k != "geometry"},
        "source_crs": source.crs.to_string(),
        "source_crs_wkt": source.crs.to_wkt(),
        "map_crs": CRS.from_user_input(target_crs).to_string(),
        "transformation": "GeoPandas/pyproj to_crs; explicit CRS, x=easting, y=northing",
        "pyproj_operation": operation.description,
        "pyproj_pipeline": operation.definition,
        "pyproj_operation_accuracy_m": operation.accuracy,
        "original_bounds": geometry.bounds,
        "map_bounds_m": transformed.geometry.iloc[0].bounds,
        "map_area_m2": transformed.geometry.iloc[0].area,
        "source_vertex_count": sum(len(p.exterior.coords) + sum(len(r.coords) for r in p.interiors)
                                   for p in (geometry.geoms if geometry.geom_type == "MultiPolygon" else [geometry])),
    }
    return transformed.geometry.iloc[0], metadata


def configure_font(font=None):
    family = "DejaVu Sans"
    if font is not None:
        font_manager.fontManager.addfont(str(font))
        family = font_manager.FontProperties(fname=font).get_name()
    plt.rcParams.update({"font.family": family, "font.size": 8,
        "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "legend.fontsize": 7, "pdf.fonttype": 42,
        "ps.fonttype": 42, "axes.linewidth": .6, "lines.linewidth": 1.1,
        "savefig.facecolor": "white"})


def make_figure(table, flights, boundary=None):
    """Retain the reference panel structure, colors, samples and counts."""
    recurrent = recurrence_mask(table)
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 5.6), gridspec_kw={"width_ratios": [1, 1]})
    fig.subplots_adjust(left=.10, right=.98, top=.82, bottom=.095, wspace=.55)
    ax = axs[0]
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#6F4E37", "#56B4E9"]
    markers = ["o", "s", "^", "D", "v", "P"]
    handles = []
    # Draw the actual vector boundary beneath every unchanged point/ring.
    if boundary is not None:
        polygons = boundary.geoms if boundary.geom_type == "MultiPolygon" else [boundary]
        for polygon in polygons:
            for ring in [polygon.exterior, *polygon.interiors]:
                xy = np.asarray(ring.coords)
                ax.plot(xy[:, 0] / 1000, xy[:, 1] / 1000, color=".3", lw=.7, zorder=0)
    for (date, group), color, marker in zip(table.groupby("date", sort=True), colors, markers):
        ax.scatter(group.x_m / 1000, group.y_m / 1000, s=12, color=color,
                   marker=marker, edgecolor="white", linewidth=.3)
        handles.append(Line2D([], [], ls="", marker=marker, ms=4, color=color, label=date))
    ax.scatter(table.loc[recurrent, "x_m"] / 1000,
               table.loc[recurrent, "y_m"] / 1000,
               s=38, facecolor="none", edgecolor="black", linewidth=.5)
    # Fixed marker rings identify recurrence; their size is not a 125-m buffer.
    ax.set_xlabel("UTM easting (km)")
    ax.set_ylabel("UTM northing (km)")
    ax.set_aspect("equal", adjustable="box")
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.annotate("Grid N", xy=(.9, .97), xytext=(.9, .81), xycoords="axes fraction",
                ha="center", fontsize=7,
                arrowprops={"arrowstyle": "-|>", "lw": .6, "color": ".2"})
    ax.set_title("Threshold-exceeding cells")
    ax.text(.02, .98, "(a)", transform=ax.transAxes, ha="left", va="top", fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": .8, "pad": .5})
    ax = axs[1]
    counts = table.groupby("date").size().reindex(flights, fill_value=0)
    yy = np.arange(len(flights))
    ax.barh(yy, counts.to_numpy(), color=".55", height=.55)
    ax.set_yticks(yy, flights)
    ax.set_ylim(13.7, -.7)
    ax.set_xlim(0, 44)
    ax.set_xlabel("Number of paired 50-m cells")
    ax.set_title("Counts by ASO flight")
    for y, n in zip(yy, counts):
        ax.text(n + .4, y, str(n), va="center", fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.text(.99, .98, "(b)", transform=ax.transAxes, ha="right", va="top", fontsize=8)
    handles.append(Line2D([], [], ls="", marker="o", ms=6, mfc="none", mec="black",
                          mew=.5, label="Within 125 m"))
    if boundary is not None:
        handles.append(Line2D([], [], color=".3", lw=.7, label="Watershed boundary"))
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .995),
               ncol=3, frameon=False, fontsize=7, columnspacing=1.4)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--locations", type=Path, required=True)
    parser.add_argument("--flights", type=Path, required=True,
                        help="Accepted 14-flight table with an ISO date column")
    parser.add_argument("--boundary", type=Path, required=True,
                        help="Original WBDHU10 shapefile or equivalent vector file with a declared CRS")
    parser.add_argument("--huc-field", default="huc10")
    parser.add_argument("--font", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("figures"))
    args = parser.parse_args()
    table = load_locations(args.locations)
    flights = load_flights(args.flights)
    boundary, metadata = load_boundary(args.boundary, table.crs.iloc[0], args.huc_field)
    # This is a placement audit only. Do not use it to filter or remove samples.
    covered = [boundary.covers(Point(x, y)) for x, y in table[["x_m", "y_m"]].to_numpy()]
    if not all(covered):
        raise ValueError("Reference centers and the selected watershed do not agree")
    configure_font(args.font)
    figure = make_figure(table, flights, boundary)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for ext in ("pdf", "png"):
        path = args.output_dir / f"JHM_S1_aso_extremes.{ext}"
        figure.savefig(path, dpi=600, facecolor="white")
        outputs[path.name] = sha256(path)
    plt.close(figure)
    boundary_files = [args.boundary]
    if args.boundary.suffix.lower() == ".shp":
        boundary_files = [args.boundary.with_suffix(s) for s in (".shp", ".dbf", ".shx", ".prj")]
    record = {"operation": "Overlay actual watershed vector; saved original observations unchanged; no models rerun",
        "software_versions": software_versions(),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "inputs_sha256": {str(p): sha256(p) for p in [args.locations, args.flights, *boundary_files]},
        "font_sha256": sha256(args.font) if args.font else None,
        "boundary": metadata, "record_count": len(table), "flight_count": len(flights),
        "counts_by_flight": table.groupby("date").size().reindex(flights, fill_value=0).to_dict(),
        "recurrent_records": int(recurrence_mask(table).sum()),
        "recurrence_indicator": "Fixed-size black marker rings indicate another flagged center on a different flight within 125 m (inclusive); the rings are not literal 125-m-radius buffers",
        "centers_covered_by_actual_polygon": int(sum(covered)),
        "outputs_sha256": outputs}
    (args.output_dir / "ASO_extremes_execution.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"outputs": outputs, "records": len(table), "all_centers_inside": all(covered)}, indent=2))


if __name__ == "__main__":
    main()
