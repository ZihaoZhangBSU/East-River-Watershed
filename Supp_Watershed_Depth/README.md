# Supplemental watershed snow depth

`Supp_Watershed_Depth.py` generates `JHM_S_watershed_depth.pdf/png`, the depth
counterpart of current Fig. 2. It uses the same nine water-year panels, figure
size (6.5 by 5.5 inches), colors, line styles, black ASO markers, date windows
and annotation layout. Typography matches Fig. 2: 8-point labels/titles,
7-point ticks/legend and 6.5-point panel statistics. The partial-year title uses
the ASCII spelling `WY2026 (Oct-Jun)` to avoid encoding-dependent punctuation.

Research inputs and generated figures are supplied separately from this code.
Original inputs are read only. Neither script runs a snow-model simulation.

## Build depth-specific inputs

The required existing inputs are:

- Eighteen production daily spatial caches, named
  `{isnobal,noahmp}_wyYYYY_noahgrid_daily.npz` for WY2018 through WY2026.
  Each must contain `dates` and its independently remapped `snow_depth_mm`
  array, in millimetres, on the canonical Noah-MP grid.
- `grid_common.npz`, with fractional watershed weights, fixed domain, affine
  transform and projected CRS. The accepted domain has 842 positive-area
  targets representing 748.679519619 square kilometres.
- The validated `F5_corrected_fields.npz` and original
  `F5_correction_manifest.json`. These contain the independently remapped ASO
  depth fields, their depth-specific full-target valid-area coverage and
  original ASO source hashes. SWE masks are never substituted.

Example commands use placeholders for local research input locations:

```text
python extract_watershed_depth.py --cache-dir CACHE --grid grid_common.npz --corrected-fields F5_corrected_fields.npz --corrected-manifest F5_correction_manifest.json --output-dir data
```

Optional `--cache-manifest model_operator_provenance.json` checks all eighteen
cache content hashes against the accepted source manifest. Optional
`--archived-daily core_watershed_daily.csv` checks the newly reduced depth
means against archived common-grid depth means. Use a separate output folder.

The extractor writes:

- `models/watershed_depth_daily.csv`: all 3,195 daily means and support fields.
- `models/depth_variable_coverage.csv`: depth validity/zero/negative counts.
- `aso/aso_depth_flights.csv`: all 14 weighted ASO depth means, own retained
  depth counts/areas/coverage, zero/extreme information and original hashes.
- `aso/matched_depth_flight_means.csv`: all three sources sampled on identical
  flight-specific retained depth targets, for direct matched comparisons.
- `aso/aso_depth_retained_fields.npz` and `depth_extraction_provenance.json`.

Daily model depth requires complete finite support on every fixed target;
missing depth aborts rather than changing the watershed denominator. ASO uses
its own finite depth mask and its own independently remapped full-target
coverage, retaining targets with Cj >= 0.5. The basin mean is the sum of depth
times fractional watershed area divided by retained fractional watershed area.
There is no Cj multiplier and no equal-cell weighting. Finite zeros, negative
values and extremes remain. Coverage means retained target watershed area
divided by fixed watershed area, rather than exact native flight-polygon area.

The ASO fields inherit the already verified raw-raster remapping and source
metadata. This extraction checks the frozen field archive against its original
manifest and freshly reduces depth; it does not reread/remap every original
50-m model day or original ASO raster. The production iSnobal cache remapped
depth independently of SWE with area-overlap averaging. Per-day native valid
subcell arrays are not saved, so complete target-cell depth does not establish
complete native subcell coverage. Saved model flight snapshots are required to
match production daily depth fields exactly before matched means are emitted.

## Generate the figure

```text
python Supp_Watershed_Depth.py --data-dir data --output-dir figures --font FONT.ttf
```

Omit `--font` for Matplotlib's bundled sans-serif. The accepted rendering uses
the same Nimbus Sans font file as Fig. 2. Input tables must explicitly identify
`Depth`; renaming a SWE table is rejected. The figure contains 18 daily model
curves and all 14 ASO depth markers, with no interpolation of flights.

Per-panel r, MB and relative MB compare the two daily model watershed-depth
series over exactly the dates shown. They are temporal intermodel agreement,
not ASO accuracy metrics or the per-flight spatial r used in Fig. 5. States
are at 2300 UTC. MB is mean(Noah-MP depth minus iSnobal depth), in millimetres;
relative MB is 100 times MB divided by mean iSnobal depth. Every displayed
calendar day, including valid snow-free zero, remains. Constant-series r and
zero-denominator relative MB are shown as NA rather than zero.

WY2018-WY2025 display full October-September years; WY2026 ends on 30 June.
`Watershed_depth_model_statistics.csv` contains these nine annotated periods.
`Watershed_depth_common_OctJun_statistics.csv` separately supplies common
October-June sensitivity scores; those do not replace the panel statistics.
No interannual mean-ratio panel is present. ASO flight footprints differ from
the complete fixed model domain, so direct flight comparisons should use the
matched table rather than subtracting a flight marker from a full-basin curve.

## Verify

```text
python -m unittest discover -s tests -v
```

The tests cover fractional versus equal-cell weighting, retained finite zeros
and extremes, depth-specific validity and inclusive coverage threshold,
missing fixed-domain depth, full/leap/partial calendars, wrong-variable inputs,
temporal-statistic signs/denominators, unavailable correlations and figure
completeness. Local validation additionally compared all 6,390 daily means,
nine panel-statistic rows and fourteen ASO depth means/counts/coverages with
independent saved-field reductions that import neither production script.

Dependencies: Python 3.10+, NumPy, pandas, Matplotlib and rasterio. Saved data,
raw observations and local validation receipts are not distributed here.
