# Figure 4: WY2026 SNOTEL SWE

This portable script creates the three station SWE panels at **5.85 x 5.67 inches**, with **8 pt ticks, legend, axis labels, panel letters, and station titles**. It preserves the accepted daily curves, iSnobal cell ranges, 18 timing markers, date-range shading, January–June display, and existing station asterisk convention. It replays cached tables only.

Use an existing accepted `SNOTEL_WY2026_SWE_plot_data.csv` (819 station-day rows: 273 dates per station) and `station_timing_full_record.csv`. The supplied timing dates are plotted as provided: this script does not recalculate peak or snow-disappearance definitions. A peak marker sits at its daily SWE value; a confirmed disappearance marker sits at 10 mm.

```sh
python Fig4_SNOTEL_WY2026.py \
  --plot-data /path/to/SNOTEL_WY2026_SWE_plot_data.csv \
  --timing-data /path/to/station_timing_full_record.csv \
  --font /path/to/NimbusSans-Regular.ttf \
  --output-dir outputs
```

Dependencies are listed in `requirements.txt`. The font path is explicit; inputs and font are not bundled with the public code. The script rejects incomplete calendars, duplicate dates, nonfinite/negative SWE, inconsistent range bounds, missing/duplicate WY2026 timing records, and disappearance dates on or before the corresponding source's peak. It records input hashes and software versions in `generation_receipt.json`. Run the synthetic input checks with `python -m unittest discover -s . -p 'test_fig4_inputs.py'`.

**Placement:** use `width=0.90\textwidth` for a manuscript with a 6.5-inch text block, or `width=5.85in`. Both place this PDF at its intrinsic size. Do not apply another 90% reduction to that placement. The gold band spans 18 March through the start of 22 March 2026, preserving the four inclusive dates 18–21 March.

The compact layout derives from the accepted `station()` function in the publication replay. All scientific arrays are unchanged. Only the physical canvas and the final tick/legend font sizes differ from the earlier 6.5-inch-wide figure.
