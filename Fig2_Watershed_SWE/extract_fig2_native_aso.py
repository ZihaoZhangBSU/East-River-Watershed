#!/usr/bin/env python3
"""Extract Fig. 2 SWE means on native, unfiltered ASO flight footprints.

Requires numpy, rasterio, shapely >=2, pyproj, and pyogrio. Source inputs are
opened read-only. No SWE magnitude/positivity screen or model-grid remapping is
applied. Means use exact planar source-pixel/watershed intersections in the
projected native raster CRS, including fractional boundary pixels. Coverage is
valid overlap area / FULL watershed polygon area (including areas outside the
raster extent). Areas are converted from native projected units to m2.

Example:
  python extract_fig2_native_aso.py --aso-directory ASO --watershed WBDHU10.shp \
    --watershed-field name --watershed-value "East River" --flights flights.csv \
    --output native_ASO_14_flights.csv --manifest native_ASO_manifest.json

The flights CSV needs an ISO-8601 ``date`` column. The optional ``source_file``
column overrides --filename-pattern for a flight. Source paths are interpreted
relative to --aso-directory. The default ASO file metadata must identify SWE
units as metres or millimetres; --swe-units can specify known source units.
Geographic rasters and rotated grids are rejected explicitly. Do not use square
degrees as physical area; reproject such sources with an appropriate validated
workflow first. Planimetric projected area is a map-area convention, not a
geodesic ground-area estimate.
"""

from __future__ import annotations

import argparse
import csv
from datetime import date
import hashlib
import json
from pathlib import Path
import platform
import re
import sys

sys.dont_write_bytecode = True

import numpy as np
import pyogrio
from pyogrio.raw import read as read_vector
from pyproj import CRS, Transformer
import rasterio
from rasterio.features import rasterize
from rasterio.windows import Window
import shapely
from shapely.ops import transform as transform_geometry


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def select_watershed(path: Path, field: str, value: str):
    info = pyogrio.read_info(path)
    matches = [str(name) for name in info["fields"] if str(name).casefold() == field.casefold()]
    if len(matches) != 1:
        raise ValueError(f"Expected a unique field {field!r}; available: {list(info['fields'])}")
    metadata, _, geometries, columns = read_vector(path, columns=[matches[0]])
    selected = [i for i, item in enumerate(columns[0])
                if str(item).strip().casefold() == value.strip().casefold()]
    if len(selected) != 1:
        raise ValueError(f"Expected exactly one exact watershed match for {value!r}; got {len(selected)}")
    geometry = shapely.from_wkb(geometries[selected[0]])
    if geometry is None or geometry.is_empty or not geometry.is_valid:
        raise ValueError("The selected watershed is empty or invalid; repair it explicitly before extraction")
    if geometry.geom_type not in {"Polygon", "MultiPolygon"}:
        raise ValueError("The watershed must be Polygon or MultiPolygon")
    if not metadata["crs"]:
        raise ValueError("The watershed has no declared CRS")
    return geometry, CRS.from_user_input(metadata["crs"])


def projected_area_factor(crs: CRS) -> float:
    """Return native area-to-m2 factor, refusing invalid area coordinates."""
    if not crs.is_projected:
        raise ValueError("Native ASO raster must have a projected CRS; geographic degree areas are unsupported")
    if len(crs.axis_info) < 2:
        raise ValueError("Raster CRS has no two-dimensional projected axes")
    factors = [float(axis.unit_conversion_factor) for axis in crs.axis_info[:2]]
    if not np.all(np.isfinite(factors)) or min(factors) <= 0:
        raise ValueError("Raster CRS projected-axis unit conversion is not valid")
    return factors[0] * factors[1]


def swe_units_factor(src, requested: str) -> tuple[float, str]:
    if requested == "m":
        return 1000.0, "m (explicit CLI)"
    if requested == "mm":
        return 1.0, "mm (explicit CLI)"
    evidence = []
    if src.units[0]:
        evidence.append(str(src.units[0]))
    unit_keys = {"unit", "units", "unittype", "swe_unit", "swe_units",
                 "snow_water_equivalent_unit", "snow_water_equivalent_units"}
    describes_swe = lambda text: bool(re.search(r"\bswe\b|snow[ _-]*water[ _-]*equivalent", str(text).lower()))
    if src.descriptions[0] and describes_swe(src.descriptions[0]):
        evidence.append(str(src.descriptions[0]))
    for tags in (src.tags(), src.tags(1)):
        for key, value in tags.items():
            if key.lower() in unit_keys or describes_swe(value):
                evidence.append(str(value))
    identified = set()
    for text in evidence:
        text = text.lower()
        if re.search(r"\bmm\b|\bmillimet(?:er|re)s?\b", text):
            identified.add("mm")
        if re.search(r"\bm\b|\bmet(?:er|re)s?\b", text):
            identified.add("m")
    if len(identified) > 1:
        raise ValueError(f"Conflicting SWE unit metadata for {src.name}: {evidence}")
    if not identified:
        raise ValueError(f"Cannot establish SWE units from metadata for {src.name}; supply --swe-units m or mm")
    unit = next(iter(identified))
    return (1000.0 if unit == "m" else 1.0), f"{unit} (source metadata)"


def validate_output_paths(output: Path, manifest: Path, input_paths):
    """Refuse every output/input collision before creating or overwriting files."""
    resolved_output, resolved_manifest = output.resolve(), manifest.resolve()
    if resolved_output == resolved_manifest:
        raise ValueError("CSV output and manifest must be distinct paths")
    protected = {path.resolve() for path in input_paths}
    if resolved_output in protected or resolved_manifest in protected:
        raise ValueError("Output paths must not overwrite any source input")


def overlap_weights(src, watershed, watershed_crs: CRS):
    """Return local read window, native pixel overlap areas and full area in m2."""
    crs = CRS.from_user_input(src.crs) if src.crs else None
    if crs is None:
        raise ValueError("ASO raster has no declared CRS")
    area_factor = projected_area_factor(crs)
    affine = src.transform
    if affine.b != 0 or affine.d != 0 or affine.a <= 0 or affine.e >= 0:
        raise ValueError("Only unrotated north-up ASO grids are supported")
    local = watershed if watershed_crs == crs else transform_geometry(
        Transformer.from_crs(watershed_crs, crs, always_xy=True).transform, watershed)
    full_area = float(local.area * area_factor)
    if not np.isfinite(full_area) or full_area <= 0:
        raise ValueError("Watershed area must be finite and positive")
    left, bottom, right, top = local.bounds
    c0 = max(0, int(np.floor((left - affine.c) / affine.a)))
    c1 = min(src.width, int(np.ceil((right - affine.c) / affine.a)))
    r0 = max(0, int(np.floor((top - affine.f) / affine.e)))
    r1 = min(src.height, int(np.ceil((bottom - affine.f) / affine.e)))
    if c1 <= c0 or r1 <= r0:
        raise ValueError("The source raster does not overlap the watershed")
    window = Window(c0, r0, c1 - c0, r1 - r0)
    transform = src.window_transform(window)
    touched = rasterize([(shapely.geometry.mapping(local), 1)],
                        out_shape=(r1 - r0, c1 - c0), transform=transform,
                        all_touched=True, fill=0, dtype="uint8")
    rows, cols = np.nonzero(touched)
    weights = np.zeros(touched.shape, dtype="float64")
    cell_area = abs(affine.a * affine.e)
    shapely.prepare(local)  # Accelerates repeated covers predicates; geometry is unchanged.
    for start in range(0, len(rows), 50000):
        row, col = rows[start:start + 50000], cols[start:start + 50000]
        x0 = transform.c + col * affine.a
        y0 = transform.f + row * affine.e
        boxes = shapely.box(x0, y0 + affine.e, x0 + affine.a, y0)
        inside = shapely.covers(local, boxes)
        areas = np.full(len(row), cell_area, dtype="float64")
        partial = ~inside
        if np.any(partial):
            areas[partial] = shapely.area(shapely.intersection(boxes[partial], local))
        weights[row, col] = np.clip(areas, 0.0, cell_area) * area_factor
    return window, weights, full_area


def calculate_mean(values: np.ndarray, overlap_m2: np.ndarray, full_area_m2: float):
    """SWE may be zero, negative or very large: exclude only missing/nonfinite."""
    valid = np.isfinite(values) & (overlap_m2 > 0)
    valid_area = float(np.sum(overlap_m2[valid], dtype="float64"))
    if valid_area <= 0:
        raise ValueError("No finite native ASO SWE intersects the watershed")
    fraction = valid_area / full_area_m2
    if fraction > 1 + 1e-9:
        raise ValueError(f"Valid source area exceeds full watershed area: {fraction}")
    return float(np.sum(values[valid] * overlap_m2[valid], dtype="float64") / valid_area), min(1.0, fraction), valid_area


def extract(aso_directory: Path, watershed_path: Path, flights_path: Path,
            output_path: Path, manifest_path: Path, watershed_field="name",
            watershed_value="East River", swe_units="auto",
            filename_pattern="ASO_50M_SWE_USCOGE_{yyyymmdd}.tif"):
    watershed, watershed_crs = select_watershed(watershed_path, watershed_field, watershed_value)
    with flights_path.open(newline="", encoding="utf-8-sig") as stream:
        flights = list(csv.DictReader(stream))
    dates = [date.fromisoformat(item["date"]) for item in flights]
    if len(set(dates)) != len(dates) or not dates:
        raise ValueError("Flight dates must be nonempty and unique")
    ordered_flights = sorted(zip(dates, flights), key=lambda pair: pair[0])
    source_files = [aso_directory / (flight.get("source_file") or filename_pattern.format(
        date=day.isoformat(), yyyymmdd=day.strftime("%Y%m%d"))) for day, flight in ordered_flights]
    siblings = sorted(watershed_path.parent.glob(watershed_path.stem + ".*")) if watershed_path.suffix.lower() == ".shp" else [watershed_path]
    validate_output_paths(output_path, manifest_path, [flights_path, watershed_path, *siblings, *source_files, Path(__file__)])
    cache, rows, sources = {}, [], []
    for (day, flight), source in zip(ordered_flights, source_files):
        source_hash = sha256(source)
        with rasterio.open(source) as src:
            if src.count != 1:
                raise ValueError(f"Expected single-band SWE raster: {source}")
            key = (src.crs.to_wkt() if src.crs else "", tuple(src.transform), src.width, src.height)
            if key not in cache:
                cache[key] = overlap_weights(src, watershed, watershed_crs)
            window, weights, full_area = cache[key]
            multiplier, unit_evidence = swe_units_factor(src, swe_units)
            # The raster's validity mask encodes NoData/explicit source masks.
            # GDAL scale and offset are separate from a masked read.
            values = src.read(1, window=window, masked=True).astype("float64").filled(np.nan)
            values = (values * src.scales[0] + src.offsets[0]) * multiplier
            mean, fraction, valid_area = calculate_mean(values, weights, full_area)
            source_details = {
                "date": day.isoformat(), "source_file": source.name, "sha256": source_hash,
                "crs": src.crs.to_string(), "crs_wkt": src.crs.to_wkt(),
                "shape": [src.height, src.width], "transform": list(src.transform),
                "nodata": (str(src.nodata) if src.nodata is not None and not np.isfinite(src.nodata) else src.nodata),
                "scale": src.scales[0], "offset": src.offsets[0],
                "mask_flags": [[flag.name for flag in flags] for flags in src.mask_flag_enums],
                "units": unit_evidence, "source_tags": src.tags(),
                "covered_grid_area_m2": float(weights.sum()),
                "raster_extent_watershed_area_fraction": float(weights.sum() / full_area),
            }
        rows.append({"date": day.isoformat(), "wy": day.year + (day.month >= 10),
                     "aso_watershed_mean_mm": mean, "aso_valid_area_fraction": fraction,
                     "valid_area_m2": valid_area, "watershed_area_m2": full_area,
                     "source_file": source.name, "source_sha256": source_hash})
        sources.append(source_details)
        print(f"{day}: {mean:.9f} mm; full-watershed valid coverage {100*fraction:.6f}%", flush=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "algorithm": "native ASO exact planar pixel-watershed intersections; no magnitude/positivity screen; finite source-valid SWE only",
        "coverage_denominator": "Full projected watershed polygon area, including areas beyond each raster extent",
        "area_convention": "Native projected-CRS planimetric area, converted to square metres; not geodesic ground area",
        "mean_equation": "sum(finite_SWE_mm * valid_pixel_watershed_overlap_m2) / sum(valid_pixel_watershed_overlap_m2)",
        "watershed_field": watershed_field, "watershed_value": watershed_value,
        "watershed_crs": watershed_crs.to_string(),
        "watershed_sources": [{"file": item.name, "sha256": sha256(item)} for item in siblings],
        "flight_manifest": {"file": flights_path.name, "sha256": sha256(flights_path)},
        "output": {"file": output_path.name, "sha256": sha256(output_path)},
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "rasterio": rasterio.__version__, "shapely": shapely.__version__,
                     "pyogrio": pyogrio.__version__, "pyproj": __import__("pyproj").__version__},
        "sources": sources,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return rows, manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--aso-directory", type=Path, required=True)
    parser.add_argument("--watershed", type=Path, required=True)
    parser.add_argument("--watershed-field", default="name")
    parser.add_argument("--watershed-value", default="East River")
    parser.add_argument("--flights", type=Path, required=True)
    parser.add_argument("--filename-pattern", default="ASO_50M_SWE_USCOGE_{yyyymmdd}.tif")
    parser.add_argument("--swe-units", choices=["auto", "m", "mm"], default="auto")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    extract(args.aso_directory, args.watershed, args.flights, args.output, args.manifest,
            args.watershed_field, args.watershed_value, args.swe_units, args.filename_pattern)


if __name__ == "__main__":
    main()
