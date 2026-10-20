# Station precipitation/SWE bias and storage/density diagnostics

`Station_Precip_SWE_Bias_and_Storage_Density.py` produces two supplemental figures:

- `JHM_S_station_precip_swe_bias.pdf` and `.png`: rows are precipitation mean bias (mm/day), precipitation relative mean bias (%), SWE mean bias (mm), and SWE relative mean bias (%); columns are Butte, Schofield Pass, and Upper Taylor*. Models are iSnobal (blue circles) and Noah-MP (orange dashed squares).
- `JHM_S_station_storage_density.pdf` and `.png`: rows are own-peak SWE / own-source precipitation and annual mean daily density; columns use the same stations. SNOTEL is black, iSnobal blue, and Noah-MP orange. Shading marks partial WY2026.

These reproduce the diagnostic types in the older manuscript's S2 and S3 using the accepted current calculation periods. The current supplement's existing S1-S11 must remain in place; the bias and storage/density figures append as S12 and S13. Both assets are 6.5 inches wide, with 8 pt ticks and labels. Operator/period/gap-fill explanations belong in the captions. The Upper Taylor* annotation is retained; grey shading identifies partial WY2026.

## Run

Install `requirements.txt`, then run:

```text
python Station_Precip_SWE_Bias_and_Storage_Density.py --data-dir INPUT_TABLES --output-dir OUTPUTS --font NimbusSans-Regular.ttf --fig3-score-csv SNOTEL_figure_data.csv
python -m unittest -v test_station_diagnostics
```

The `--font` and `--fig3-score-csv` arguments are optional. The supplied journal renders use Nimbus Sans. Without that font, the program uses DejaVu Sans. No absolute private paths are embedded in the code. The journal font is supplied separately and is not distributed here. The Fig. 3 cross-check compares all 54 SWE biases, paired day counts, and period endpoints against the accepted score table; the figure generation fails if they differ.

`INPUT_TABLES` must contain the original frozen `snotel_daily.csv`, `isnobal_station_daily.csv`, `noahmp_station_daily.csv`, and `snotel_precip_gapfill_log.csv`. Research inputs, generated PDFs/PNGs, and diagnostic CSVs are not part of the public code upload. The script records input, font, code, and output SHA256 hashes and runtime versions in `generation_receipt.json`. No model is rerun and no input is edited.

## Definitions and validation

All three sources must include exactly one value on every expected station/date: full October-September WY2018-WY2025 and October-June WY2026. Leap days and zero values remain. Input validation rejects duplicate/extra/missing dates, inconsistent year labels, missing required SWE/depth/precipitation values, or negative SWE/depth/precipitation. These archived extracts have 3,195 dates per station. Station observations use `swe_mm`, `snow_depth_mm`, `precip_gapfilled_mm`; iSnobal uses the corresponding `*_nearest_mm` fields; Noah-MP uses `swe_mm`, `snow_depth_mm`, `precip_mm`.

Bias is mean(model - SNOTEL) over every paired day. Relative bias is `100 * mean(model - SNOTEL) / mean(SNOTEL)`, not a mean of daily percentage errors. Precipitation values are daily increments, so precipitation mean bias has units mm/day; the annual total precipitation difference equals this mean bias times the paired-day count. SWE mean bias is in mm. An all-zero reference mean produces an explicit undefined relative-bias status. This case does not occur in the real inputs.

Each source's peak is the first exact maximum of its unsmoothed SWE between 1 October and 30 June. The denominator is that same source's complete daily precipitation sum from 1 October through its own peak date, inclusive. Missing precipitation or a missing calendar date is an error. A zero precipitation denominator produces an explicit undefined ratio status. An annual SWE maximum after June does not change the accepted October-June peak anchor.

Density is the arithmetic mean of daily `1000 * SWE_mm / depth_mm` (kg/m3), using finite SWE/depth pairs with depth **at least 100 mm**, separately for each source. It uses the entire available annual period stated above. There is no common source/date mask and no density upper cap. Source-specific accepted day counts are exported. No valid support yields an explicit undefined status. Consequently the density means can reflect different snow-season dates; they are not a matched-date density-error score.

SNOTEL precipitation uses four archived replacements with nearest iSnobal precipitation: Butte 2025-06-07 and 2025-06-13; Upper Taylor 2021-10-02 and 2021-12-17 (WY2022). The program verifies the daily missing-observation flags, source labels, replacement values against iSnobal, and the complete gap-fill log. It performs no new imputation. Model precipitation is never replaced. The filled reference is not fully independent of iSnobal. Annual replacement counts/sums/dates and replacement counts through each SNOTEL peak are exported.

The storage/precipitation ratio combines source-specific peak dates and windows; it is not a snow-retention efficiency. Annual precipitation amount, timing, phase, snow redistribution, and energy processes cannot be uniquely separated by these panels. Positive SWE bias with small or negative precipitation bias motivates further diagnosis; it does not identify a unique cause.

## Numeric outputs

- `station_annual_precip_swe_bias.csv`: 108 station/year/model/variable records with biases, paired dates, reference/model means and sums, and gap-fill support.
- `station_annual_storage_density.csv`: 81 station/year/source records with peak dates, peak SWE, own-peak precipitation sum/day count, ratios, density means and accepted/excluded day counts.
- `station_daily_diagnostic_support.csv`: 28,755 source/station/dates, retaining all input values, source-specific density flags and values, and SNOTEL-only precipitation replacement flags.
- `station_precip_gapfill_support.csv`: 27 station/year records with replacement counts/sums/dates.
- `station_precip_gapfill_dates.csv`: four original replacement records.

The synthetic tests cover zeros and denominator definitions, leap-day and partial-year support, incomplete precipitation/date failures, tied/own-source peak dates, post-June SWE maxima, the inclusive 100 mm density threshold, daily-mean versus ratio-of-means density, explicit no-support cases, absence of a density cap, and replacement support through the peak.
