Fig. 5 (SWE) and supplemental ASO depth spatial performance

Fig5_ASO_Performance.py recomputes all fourteen SWE/depth flight comparisons and
the lower/middle/upper elevation summaries. It shares the portable validated
extract_fig2_common_aso.py remapping module at the figure folder's top level.
The plotter emits 6.5-inch-wide figures with 8-point text:
  JHM_F5_aso_performance.pdf/png: SWE MB (a), RMSE (b), spatial Pearson r (c).
  JHM_S_ASO_depth_performance.pdf/png: depth MB (a) and RMSE (b).
Both use the same fourteen flights, unchanged full-precision scores and paired
supports, colors, ASO-mean marker areas, and WY2026 shading. MB is model minus
ASO. The split changes presentation only; it does not change a score or run a
snow-model simulation. The depth figure's supplement number is set by the paper.

Replay already computed, corrected full-precision scores without reading raw
rasters (complete canonical F5_elevation_metrics_corrected.csv, 224 rows):
  python code/f5/Fig5_ASO_Performance.py \
    --scores-csv data/f5_boundary/F5_elevation_metrics_corrected.csv \
    --output-directory data/f5_replay --figure-directory figures \
    --font assets/NimbusSans-Regular.ttf --figure-selection swe
Use --figure-selection swe for a main-figure-only revision; it never writes
the accepted supplemental depth PDF/PNG. The default both retains the earlier
two-figure workflow, and depth selects only the unchanged two-panel depth plot.
The depth plot's numerical values and artist geometry have not changed.
Replay validates the exact schema, all fourteen flight dates, the complete
model/variable/elevation keys, finite required scores, relative-score consistency,
and elevation support-count partitions. It also refuses common marker sizing
for unequal model reference supports. Replay validation checks the table's
internal consistency; it does not establish raw-raster provenance. Preserve the
original corrected generator manifest and independently verify against saved
fields with verify_f5_numeric.py. The replay manifest records the input hash,
generator hash, runtime versions, optional font hash, and final asset hashes.
Raw-raster inputs and --baseline-displayed-csv cannot be mixed with --scores-csv.

Required local inputs (not distributed):
  * Original ASO_50M_{SWE,SD}_USCOGE_YYYYMMDD.tif files for fourteen flights.
  * flights.csv with date column (ISO dates).
  * Canonical grid_common.npz with weights, fixed_domain_mask, transform, crs,
    x/y and finite target DEM (dem).
  * Saved unfiltered D23 model rasters on the Noah grid:
    YYYYMMDD_{swe,snow_depth}_{isnobal_on_noah,noah_native}.tif.
No new snow-model simulations are performed. Inputs are opened read-only.

Example (paths are placeholders for the user's actual local data):
  python code/f5/Fig5_ASO_Performance.py --aso-directory ASO \
    --model-directory model_rasters --flights flights.csv \
    --grid-weights grid_common.npz --output-directory data/f5_boundary \
    --figure-directory figures --font assets/NimbusSans-Regular.ttf

Omit --font to use Matplotlib's bundled DejaVu Sans. An optional
--baseline-displayed-csv compares the new raw metrics with a previous displayed
metrics CSV (date,variable,model,aso_mean_mm,mb_mm,rmse_mm).

Dependencies: Python 3.10+, numpy, rasterio, affine, matplotlib.
Run synthetic scoring tests:
  python -m unittest discover -s code/f5/tests -v
Run shared remapping tests at the figure folder's tests/test_extract_fig2_common_aso.py.
Run independent numeric and plot-input verification after generation:
  python code/f5/verify_f5_numeric.py --data-directory data/f5_boundary
This independently checks all 224 model/variable/elevation rows, including
paired counts, reference means, MB, RMSE, Pearson r, RMB, RRMSE, and all selected
plotted marker positions (84 for SWE-only; 140 for both figures).

Operator and interpretation
  SWE/depth honor their own NoData/masks/nonfinite values; finite zeros,
  negative and extreme values remain. ASO is independently averaged onto Noah
  targets over valid contributing source area. Its 0/1 validity raster is
  zero-padded across the full target envelope before averaging, so C_j>=0.5
  uses the FULL target footprint even at source raster edges.
  F5 retains equal-cell, model-specific finite ASO--model pair scores. Watershed
  fractions define target eligibility (w_j>0), not score weights. These scores
  differ from Fig. 2's fractional-area basin means. No C_j multiplier is used.
  The generator refuses silently assigning one ASO mean/marker size to models
  whose paired reference supports differ; the actual fourteen model pairs have
  coincident support. Saved model states are D23 UTC, not exact ASO flight time.
  Spatial r is Pearson correlation across the paired target-cell ASO and model
  SWE values separately for each flight and model. It is neither a daily
  temporal correlation nor r squared. No flights are pooled. An undefined
  correlation (constant field or insufficient pairs) is recorded as missing
  and marked NA, never zero. Exactly zero is a valid plotted correlation.
  The r axis displays the finite score range with padding bounded by [-1,1],
  with explicit labelled ticks. SWE MB/RMSE values, areas, colors and supports
  are unchanged by adding the third panel.

Outputs include two vector PDFs/600-dpi PNGs, full precision all-flight/elevation
scores, coverage audit, legacy-unpadded versus corrected metric differences,
optional previous-displayed-value comparisons, source/output hashes and figure
marker placement records. Preserve the original figure/table assets separately.
