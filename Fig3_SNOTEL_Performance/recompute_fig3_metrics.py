#!/usr/bin/env python3
"""Recompute Fig. 3 station metrics using each year's available amount record.

Historical amount scores use October--September; WY2026 uses October--June.
Peak anchors remain October--June, and SDD uses the full available record.
This script validates the saved station series; it does not extract or validate
raw observations, quality flags, raw model files, spatial lookup, or clocks.
Only Python's standard library is required.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import date, timedelta
from pathlib import Path

STATIONS = {380: "Butte", 737: "Schofield Pass", 1141: "Upper Taylor"}
YEARS = tuple(range(2018, 2027))
SOURCE_COLUMNS = {
    "SNOTEL": "swe_mm",
    "iSnobal": "swe_nearest_mm",
    "Noah-MP": "swe_mm",
}
METRIC_FIELDS = (
    "station_id", "station_name", "wy", "model", "mb", "rmse",
    "peak_date", "sdd_date", "peak_error_days", "sdd_error_days",
    "sdd_censored", "available_record_end", "observed_sdd_censored",
    "amount_window_start", "amount_window_end", "amount_n_days",
)


def calendar_days(start: date, end: date):
    """Yield every calendar date including both endpoints."""
    while start <= end:
        yield start
        start += timedelta(days=1)


def available_end(wy: int) -> date:
    return date(wy, 6, 30) if wy == 2026 else date(wy, 9, 30)


def amount_days(wy: int) -> list[date]:
    """Full historical water year; actual October--June record for WY2026."""
    return list(calendar_days(date(wy - 1, 10, 1), available_end(wy)))


def load_source(path: Path, source: str) -> dict:
    """Validate the complete canonical Fig. 3 daily source and return SWE maps.

    The strict input contract is three stations, WY2018--2026, complete daily
    calendars (historical October--September, WY2026 October--June), finite
    nonnegative SWE in mm, and unique station/date keys. Extra fields in the
    archived daily tables are accepted but do not enter these metrics.
    """
    column = SOURCE_COLUMNS[source]
    groups = {}
    required = {"date", "water_year", "station_id", "station_name", column}
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        headers = reader.fieldnames or []
        if len(set(headers)) != len(headers) or not required.issubset(headers):
            raise ValueError(f"{source}: missing or duplicate required CSV headers")
        for line, row in enumerate(reader, 2):
            where = f"{source}, CSV row {line}"
            if None in row or any(row[k] is None for k in required):
                raise ValueError(f"{where}: malformed CSV row")
            try:
                station = int(row["station_id"])
                wy = int(row["water_year"])
                day = date.fromisoformat(row["date"])
                swe = float(row[column])
            except (TypeError, ValueError) as error:
                raise ValueError(f"{where}: invalid station, year, date, or SWE") from error
            if row["date"] != day.isoformat():
                raise ValueError(f"{where}: dates must use YYYY-MM-DD")
            if station not in STATIONS or wy not in YEARS:
                raise ValueError(f"{where}: unsupported station or water year")
            if row["station_name"] != STATIONS[station]:
                raise ValueError(f"{where}: station ID/name mismatch")
            start, end = date(wy - 1, 10, 1), available_end(wy)
            if not start <= day <= end:
                raise ValueError(f"{where}: date outside the supported WY calendar")
            if not math.isfinite(swe) or swe < 0:
                raise ValueError(f"{where}: SWE must be finite and nonnegative")
            daily = groups.setdefault((station, wy), {})
            if day in daily:
                raise ValueError(f"{where}: duplicate station/date")
            daily[day] = swe
    expected_groups = {(station, wy) for station in STATIONS for wy in YEARS}
    if set(groups) != expected_groups:
        raise ValueError(f"{source}: expected all 27 station/water-year groups")
    for (station, wy), daily in groups.items():
        expected = set(calendar_days(date(wy - 1, 10, 1), available_end(wy)))
        if set(daily) != expected:
            raise ValueError(f"{source}: incomplete calendar for station {station}, WY{wy}")
    return groups


def first_peak(daily: dict, wy: int) -> tuple[date, float, list[date]]:
    """First exact maximum in October--June, retaining zero-SWE days."""
    season = {
        day: value for day, value in daily.items()
        if date(wy - 1, 10, 1) <= day <= date(wy, 6, 30)
        and value is not None and math.isfinite(value)
    }
    if not season:
        raise ValueError("No finite October--June SWE values for peak detection")
    maximum = max(season.values())
    ties = sorted(day for day, value in season.items() if value == maximum)
    return ties[0], maximum, ties


def snow_disappearance(daily: dict, peak: date) -> date | None:
    """First post-own-peak start of seven calendar days with SWE strictly <10.

    Missing dates and nonfinite values break confirmation. The available record
    end is never substituted for a confirmed event. A peak day cannot be SDD.
    """
    if not daily:
        return None
    end = max(daily)
    for start in calendar_days(peak + timedelta(days=1), end - timedelta(days=6)):
        values = [daily.get(start + timedelta(days=offset)) for offset in range(7)]
        if all(value is not None and math.isfinite(value) and value < 10.0
               for value in values):
            return start
    return None


def amount_metrics(observed: list[float], modeled: list[float]) -> tuple[float, float]:
    """Mean(model - observation) and RMSE on identical finite paired dates."""
    if len(observed) != len(modeled) or not observed:
        raise ValueError("Amount metrics require equal nonempty paired series")
    if not all(math.isfinite(value) for value in observed + modeled):
        raise ValueError("Amount metrics require finite paired values")
    differences = [model - obs for obs, model in zip(observed, modeled)]
    bias = math.fsum(differences) / len(differences)
    rmse = math.sqrt(math.fsum(value * value for value in differences) / len(differences))
    return bias, rmse


def recompute(sources: dict) -> tuple[list[dict], dict]:
    """Build the 54 plotted rows and explicit event/support audit."""
    rows, events, support = [], [], []
    for station, name in STATIONS.items():
        for wy in YEARS:
            key = (station, wy)
            days = amount_days(wy)
            daily = {source: sources[source][key] for source in SOURCE_COLUMNS}
            # All three sources must have exactly the same complete scoring
            # dates. No model-specific filtering or observation gap filling.
            for source, values in daily.items():
                if not all(day in values and math.isfinite(values[day]) for day in days):
                    raise ValueError(f"{source}: incomplete finite paired scoring dates")
            peaks = {source: first_peak(values, wy) for source, values in daily.items()}
            sdds = {source: snow_disappearance(values, peaks[source][0])
                    for source, values in daily.items()}
            support.append({"station_id": station, "wy": wy,
                            "paired_scoring_days": len(days),
                            "amount_window_start": days[0].isoformat(),
                            "amount_window_end": days[-1].isoformat(),
                            "amount_window": "october_june_partial" if wy == 2026 else "october_september_full",
                            "observed_zero_swe_days": sum(daily["SNOTEL"][d] == 0 for d in days)})
            for source in SOURCE_COLUMNS:
                peak, maximum, ties = peaks[source]
                sdd = sdds[source]
                events.append({
                    "station_id": station, "station_name": name, "wy": wy,
                    "source": source, "peak_date": peak.isoformat(),
                    "peak_swe_mm": maximum, "peak_tie_count": len(ties),
                    "peak_last_tie_date": ties[-1].isoformat(),
                    "sdd_date": sdd.isoformat() if sdd else "",
                    "sdd_confirmation_date": (sdd + timedelta(days=6)).isoformat() if sdd else "",
                    "sdd_censored": sdd is None,
                    "available_record_end": max(daily[source]).isoformat(),
                    "observed_sdd_censored": sdds["SNOTEL"] is None,
                    "amount_window_start": days[0].isoformat(),
                    "amount_window_end": days[-1].isoformat(),
                    "amount_n_days": len(days),
                })
                if source == "SNOTEL":
                    continue
                mb, rmse = amount_metrics([daily["SNOTEL"][d] for d in days],
                                          [daily[source][d] for d in days])
                valid_sdd_error = sdd is not None and sdds["SNOTEL"] is not None
                rows.append({
                    "station_id": station, "station_name": name, "wy": wy,
                    "model": source, "mb": mb, "rmse": rmse,
                    "peak_date": peak.isoformat(), "sdd_date": sdd.isoformat() if sdd else "",
                    "peak_error_days": (peak - peaks["SNOTEL"][0]).days,
                    "sdd_error_days": (sdd - sdds["SNOTEL"]).days if valid_sdd_error else "",
                    "sdd_censored": sdd is None,
                    "available_record_end": max(daily[source]).isoformat(),
                    "observed_sdd_censored": sdds["SNOTEL"] is None,
                    "amount_window_start": days[0].isoformat(),
                    "amount_window_end": days[-1].isoformat(),
                    "amount_n_days": len(days),
                })
    return rows, {"paired_support": support, "source_events": events}


def compare_reference(rows: list[dict], path: Path) -> dict:
    """Check every plotted metric and date against an optional frozen table."""
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if not set(METRIC_FIELDS).issubset(reader.fieldnames or []):
            raise ValueError("Reference must include revised amount-score window/count and censor metadata")
        reference = list(reader)
    key_fields = ("station_id", "wy", "model")
    def key(row):
        return tuple(str(row[field]) for field in key_fields)
    lookup = {key(row): row for row in reference}
    if len(lookup) != len(reference) or set(lookup) != {key(row) for row in rows}:
        raise ValueError("Reference must contain the same 54 unique station/year/model rows")
    differences, largest = [], {field: 0.0 for field in ("mb", "rmse", "peak_error_days", "sdd_error_days")}
    for row in rows:
        old = lookup[key(row)]
        for field in largest:
            a, b = row[field], old[field]
            if str(a) == "" or str(b) == "":
                if str(a) != str(b):
                    differences.append({"key": key(row), "field": field, "computed": a, "reference": b})
                continue
            difference = abs(float(a) - float(b))
            largest[field] = max(largest[field], difference)
            tolerance = 1e-10 if field in ("mb", "rmse") else 0.0
            if not math.isfinite(float(b)) or difference > tolerance:
                differences.append({"key": key(row), "field": field, "computed": a, "reference": b})
        for field in ("station_name", "peak_date", "sdd_date", "sdd_censored", "available_record_end",
                      "observed_sdd_censored", "amount_window_start", "amount_window_end", "amount_n_days"):
            if str(row[field]) != old[field]:
                differences.append({"key": key(row), "field": field, "computed": row[field], "reference": old[field]})
    return {"metric_values_checked": 4 * len(rows), "maximum_absolute_differences": largest,
            "differences": differences, "passed": not differences,
            "reference_sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()}


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snotel", type=Path, required=True, help="daily observations, SWE in swe_mm")
    parser.add_argument("--isnobal", type=Path, required=True, help="daily nearest-cell iSnobal, SWE in swe_nearest_mm")
    parser.add_argument("--noahmp", type=Path, required=True, help="daily nearest-cell Noah-MP, SWE in swe_mm")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--reference", type=Path, help="optional frozen SNOTEL_figure_data.csv comparison")
    args = parser.parse_args(argv)
    paths = {"SNOTEL": args.snotel, "iSnobal": args.isnobal, "Noah-MP": args.noahmp}
    sources = {source: load_source(path, source) for source, path in paths.items()}
    rows, audit = recompute(sources)
    audit["code_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    audit["input_files"] = [
        {"source": source, "filename": path.name,
         "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
         "swe_column": SOURCE_COLUMNS[source],
         "rows": sum(len(daily) for daily in sources[source].values())}
        for source, path in paths.items()
    ]
    audit["definitions"] = {
        "amount": "Mean(model minus SNOTEL) and RMSE, identical finite dates in all three sources: historical October-September, WY2026 October-June; zeros and leap days retained; no gap filling or July-September 2026 padding.",
        "peak": "Own first chronological exact maximum, October-June.",
        "sdd": "First day strictly after own peak starting seven consecutive calendar days SWE <10 mm; full available WY record (historical September 30, WY2026 June 30).",
        "censoring": "An unconfirmed event is blank; SDD error is blank unless both observed and model events are confirmed.",
        "errors": "Model date minus SNOTEL date in calendar days.",
        "scope": "Saved native nearest-cell daily inputs; original dates and station names retained. No raw-data extraction, observation QC, coordinate, or measurement-clock validation.",
    }
    if args.reference:
        audit["reference_comparison"] = compare_reference(rows, args.reference)
        if not audit["reference_comparison"]["passed"]:
            raise ValueError("Recomputed metrics differ from the supplied reference table")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    destinations = [args.output_dir / "SNOTEL_figure_data.csv", args.output_dir / "fig3_metrics_audit.json"]
    if any(path.resolve() in {p.resolve() for p in paths.values()} for path in destinations):
        raise ValueError("Output files must not replace source inputs")
    if args.reference and any(path.resolve() == args.reference.resolve() for path in destinations):
        raise ValueError("Output files must not replace the comparison reference")
    with destinations[0].open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=METRIC_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    destinations[1].write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Wrote {len(rows)} metric rows and {len(audit['source_events'])} source-event checks to {args.output_dir}")


if __name__ == "__main__":
    main()
