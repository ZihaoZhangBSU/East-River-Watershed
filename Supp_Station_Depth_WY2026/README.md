# Supplemental WY2026 station depth

`Supp_Station_Depth_WY2026.py` replays the accepted depth table in exactly the
compact main Figure 4 geometry: **5.85 x 3.40 inches**, with **8 pt** station
titles, panel letters, axis labels, ticks and legend. Panel order remains
Butte, Schofield Pass, Upper Taylor*. Units remain **millimetres**.

```sh
python Supp_Station_Depth_WY2026.py \
  --plot-data /path/to/SNOTEL_WY2026_depth_plot_data.csv \
  --font /path/to/NimbusSans-Regular.ttf \
  --output-dir /path/to/generated
```

Inputs and font are supplied locally and are excluded from the public source
package. The cached table has 819 station-day rows: exactly 273 dates per
station from 1 October 2025 through 30 June 2026. Required columns are `date`,
`station_id`, `SNOTEL`, `iSnobal`, `Noah-MP`, `lo`, and `hi`. The generator
requires complete, unique calendars, finite depth/range values, ordered bounds,
and a nearest iSnobal series within its range. It never fills gaps, changes a
daily value, resamples dates, reads raw model files or runs a model.

All daily values remain in the plotted artists; the view shows January-June.
The blue envelope is the daily iSnobal minimum-maximum within the corresponding
Noah-MP cell. It represents spatial variability, not a confidence interval.
Gold shading includes all four dates 18-21 March 2026. Observation and model
colors, Noah-MP dashes, the existing station asterisk and all curves/ranges
match the accepted source. The original depth panels have **no timing markers**;
the generator does not import or infer SWE peak/disappearance markers.

Outputs are `JHM_S_station_depth.pdf`, a 600 dpi PNG, a native-size 96 dpi
preview and `station_depth_generation_receipt.json`. The receipt records
input/code/output hashes, calendar support, units, dimensions and versions.
Place the PDF at `width=5.85in` to preserve the exact size and 8 pt type.
For a 6.5-inch text block, `width=0.90\textwidth` is equivalent. Do not shrink
it again. Numbering is controlled by the manuscript, not the generator.

Run the synthetic support/artist checks with:

```sh
python -m unittest discover -s . -p test_station_depth_inputs.py
```

Dependencies are Python 3.10+, NumPy, pandas and Matplotlib, listed in
`requirements.txt`. The script is standalone and does not modify main Figure 4.
