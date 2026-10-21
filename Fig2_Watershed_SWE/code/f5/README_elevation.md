# Supplemental ASO SWE elevation performance

`Supp_ASO_SWE_Elevation.py` plots an existing, validated 224-row score CSV.
It selects exactly the 84 SWE elevation-band rows and does not reopen rasters,
pool flights, recompute metrics, or run models. Keep it beside
`Fig5_ASO_Performance.py` in `code/f5/`; the helper uses that module's strict
score-schema validator and canonical fourteen flight dates. Dependencies are
the same as that figure generator (NumPy, Matplotlib, rasterio, affine and the
shared remapping module).

```powershell
python code/f5/Supp_ASO_SWE_Elevation.py `
  --scores-csv /path/to/F5_elevation_metrics_corrected.csv `
  --output-directory /path/to/elevation-output `
  --font /path/to/NimbusSans-Regular.ttf
```

The font argument is optional; without it the helper uses Matplotlib's bundled
DejaVu Sans. Font and input files are not included in the code repository.
Absolute and relative font paths are accepted, including parent-directory
components. The resolved font registry entry must identify the supplied file.

The PDF and 600 dpi PNG have native dimensions **6.5 × 7.2 inches**, with all
visible text **8.05 pt**. Placement at the manuscript template's 6.476-inch
width gives approximately 8.02 pt text and 7.17-inch height. Preserve a width
of at least 6.460 inches to retain text of at least 8 pt. An additional 96 dpi
PNG previews the native physical size.

The six panels are iSnobal (a–c) and Noah-MP (d–f), with MB in millimetres,
RMSE in millimetres, and spatial Pearson r across columns. Each panel has
fourteen chronological flight-date rows and Lower, Middle, Upper columns.
The unchanged terrain bands are Lower <2896 m, Middle 2896–<3353 m, and
Upper ≥3353 m on the canonical Noah grid. MB is model minus ASO. Equal-cell,
model-specific finite ASO–model pairs retain valid zeros and extremes; ASO
coverage and positive watershed eligibility are unchanged from Fig. 5.

Both models share fixed limits for each metric: MB −300 to +300 mm, RMSE
0–400 mm, and r −1 to +1. A future score outside these limits causes a
validation error rather than silent colour clipping. Undefined r is gray
hatched and marked NA; finite zero is retained as a numerical value and has
the ordinary zero colour. The two undefined r cells in the verified inputs
are Noah-MP Lower on 24 May 2018 and 10 June 2019, where all 187 model values
are zero. MB and RMSE remain defined for these cells.

Outputs include the exact original scalar strings and support counts in an
84-row CSV and a JSON manifest containing full plotted arrays, NA masks,
shared scales, input/code/font hashes, software versions and output hashes.
The manifest records provenance; it does not independently validate the
upstream raster operator. A separate numerical audit against saved fields is
required before publication.

Run the presentation safety tests with:

```powershell
python -m unittest discover -s code/f5 -p test_supp_aso_swe_elevation.py
```

Research inputs and generated figure assets are deliberately excluded from
the source-only repository.
