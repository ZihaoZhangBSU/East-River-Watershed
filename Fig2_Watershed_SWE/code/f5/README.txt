Fig. 5: unfiltered ASO spatial performance and full-target coverage correction

Fig5_ASO_Performance.py recomputes all fourteen SWE/depth flight comparisons and
the lower/middle/upper elevation summaries. It shares the portable validated
extract_fig2_common_aso.py remapping module at the figure folder's top level.

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
  python code/f5/tests/test_f5.py
Run shared remapping tests at the figure folder's tests/test_extract_fig2_common_aso.py.
Run independent numeric and plot-input verification after generation:
  python code/f5/verify_f5_numeric.py --data-directory data/f5_boundary
This independently checks all 224 model/variable/elevation rows, including
paired counts, reference means, MB, RMSE, Pearson r, RMB, RRMSE, and all 112
plotted marker positions.

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

Outputs include vector PDF/600-dpi PNG, full precision all-flight/elevation
scores, coverage audit, legacy-unpadded versus corrected metric differences,
optional previous-displayed-value comparisons, source/output hashes and figure
marker placement records. Preserve the original figure/table assets separately.
