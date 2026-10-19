Fig. 2 model operator and reproducibility

Run extract_fixed_model_means.py on the eighteen EXISTING production spatial
caches, not on native-grid basin means:

python extract_fixed_model_means.py --cache-dir CACHE_DIRECTORY --grid grid_common.npz --output-dir OUTPUT_DIRECTORY --native-daily figure04_native_watershed_daily.csv --archived-common-daily core_watershed_daily.csv

The last two arguments are optional comparison/check inputs. They are not used
to compute the revised model means. Required Python packages: numpy, pandas.
Run tests with: python -m unittest discover -s . -p test_fixed_models.py -v

Inputs:
  isnobal_wyYYYY_noahgrid_daily.npz and noahmp_wyYYYY_noahgrid_daily.npz,
  for YYYY=2018,...,2026; keys dates and swe_mm, with states shape [date,66,54].
  grid_common.npz: weights, fixed_domain_mask, transform, crs plus grid centers.
  Original production caches were made with specific_mass from the last daily
  iSnobal snow state and Noah-MP SNEQV at 23:00 UTC. The prior archived eighteen
  raw-date comparisons (two per year) gave zero SWE differences at the target
  cell level. The revision verifies all cached SWE cell-days and recomputes
  every basin mean, without rerunning a snow model or rereading all raw states.

Spatial method:
  iSnobal 50 m SWE states were area-overlap averaged to the native 1 km Noah-MP
  analysis grid using the verified rasterio/GDAL Resampling.average operator.
  Finite zeros are valid and no SWE magnitude/outlier filters are applied.
  Regridding averages the valid source subcells of each full target cell. A
  50 percent minimum valid-source coverage was used by the cache builder.
  Grid-edge normalization was independently checked: all 842 watershed target
  cells are entirely within the iSnobal rectangular source grid extent. This
  geometric check is distinct from native finite-data validity.

  The fixed model footprint is all target cells with positive watershed
  intersection fraction: 842 cells. This identical mask is used every date
  for both models. Each mean is sum(SWE_j*w_j)/sum(w_j), where w_j is the
  fractional watershed area of target cell j. These are not equal-cell spatial
  score weights. The weights are applied after regridding: boundary-cell SWE
  is a full-cell source average weighted by the watershed area in that cell,
  not a direct subcell integration restricted to the watershed intersection.

  The fixed-domain extractor aborts if either model contains nonfinite or
  negative SWE in any retained cell. It does not silently change the mask or
  renormalize a daily mean. Ancillary variable coverage is audited separately;
  SWE coverage is not assumed to prove mass-input or flux support.

Outputs:
  Fig2_common_grid_daily.csv: dates, water year, revised daily means and coverage.
  Fig2_old_vs_new_daily.csv / Fig2_old_vs_new_annual_metrics.csv: comparison
  against the previous native-grid means, if supplied. Annual metrics include
  each displayed period and a common October-June window separately.
  spatial_cache_variable_coverage.csv: validity, finite-zero and negative-state
  counts for each cached variable and year.
  model_operator_provenance.json: source hashes, operator, masks and checks.

Accuracy and provenance:
  Float64 accumulation is used for the supplied fractional weights and SWE.
  Independent scalar summation checks every model daily mean. Annual means and
  correlations also have independent scalar-sum checks.
  The earlier archived summary stored common valid coverage near one in
  float32 arithmetic (0.9999999403953552). The new spatial-cache check is the
  authority for exact cell-level support; the old rounded coverage alone is
  not treated as mathematical proof of complete support.

Research spatial caches, derived data, figures and manuscript text are local
handoff artifacts. The figure's public GitHub code folder contains code and
instructions only.
