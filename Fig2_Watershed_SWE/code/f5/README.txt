Fig. 5: overall SWE and depth spatial performance

Fig5_ASO_Performance.py recomputes all fourteen SWE/depth flight comparisons and
the lower/middle/upper elevation summaries. It shares the portable validated
extract_fig2_common_aso.py remapping module at the figure folder's top level.
The default plotter emits one 6.5 x 7.2-inch main figure with 8.05-point text:
  JHM_F5_aso_performance.pdf/png
  Top row: SWE MB (a), RMSE (b), spatial Pearson r (c).
  Bottom row: depth MB (d), RMSE (e), spatial Pearson r (f).
Both rows show all fourteen chronological flights and both models. All score
values, masks, colors, marker areas and WY2026 shading remain unchanged. The
new depth correlation panel uses the already validated scores. No flights are
pooled and no snow-model simulation is run. MB/RMSE units (mm) appear in panel
titles; r is dimensionless. All metric axes are shared across variable rows:
  MB [-300,200] mm, ticks -300,-150,0,150.
  RMSE [0,500] mm, ticks 0,250,500.
  Spatial r [0.65,1.0], ticks 0.7,0.8,0.9,1.0.
For future validated tables whose scores exceed these defaults, shared limits
expand to retain every finite value; valid zero/negative r are not clipped.
The main figure has separate ASO-mean marker-area legends for SWE and depth.
SWE reference=300 mm (legend samples 100,300,600); depth reference=1000 mm
(samples 300,1000,2000). Noah squares use area=16*ASOmean/reference pt^2;
open iSnobal circles use size=4*area/pi for the same enclosed marker area.
At the manuscript's 6.476-inch placement, text is at least 8.02 points.

Replay already computed, corrected full-precision scores without reading raw
rasters (complete canonical F5_elevation_metrics_corrected.csv, 224 rows):
  python code/f5/Fig5_ASO_Performance.py \
    --scores-csv data/f5_boundary/F5_elevation_metrics_corrected.csv \
    --output-directory data/f5_replay --figure-directory figures \
    --font assets/NimbusSans-Regular.ttf --figure-selection both
The default both produces only the six-panel main asset. It never writes the
former JHM_S_ASO_depth_performance.pdf/png, which should be preserved as a prior
accepted asset. The same swe/depth command-line choices remain available as
single-variable three-panel audit exports, with distinct output basenames:
  --figure-selection swe   -> JHM_F5_swe_audit.pdf/png
  --figure-selection depth -> JHM_F5_depth_audit.pdf/png
These audit-only names replace the older per-variable output names so an audit
cannot overwrite either the new main figure or the accepted depth supplement.
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
plotted marker positions (84 for a single-variable audit; 168 for the main figure
when all correlations are defined). The marker manifest counts actual nonmissing
scores, including valid zero, rather than treating unavailable r as a marker.

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
  SWE or depth values separately for each flight and model. It is neither a daily
  temporal correlation nor r squared. No flights are pooled. An undefined
  correlation (constant field or insufficient pairs) is recorded as missing
  and marked NA, never zero. Exactly zero is a valid plotted correlation.
  SWE and depth keep independent variable-specific masks and paired supports.
  Shared plot axes do not replace those masks or combine spatial observations.

Outputs include one vector PDF/600-dpi PNG, full precision all-flight/elevation
scores, coverage audit, legacy-unpadded versus corrected metric differences,
optional previous-displayed-value comparisons, source/output hashes and figure
marker placement records. Preserve the original figure/table assets separately.
