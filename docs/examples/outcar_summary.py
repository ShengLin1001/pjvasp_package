#!/usr/bin/env python3
"""Summarize selected tracked OUTCAR fixtures with the public parser classes."""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
from pathlib import Path
from tempfile import NamedTemporaryFile

from mymetal.post.newmain import PostData, PostData2, PostTime


FIELDNAMES = [
    "case",
    "status",
    "converged",
    "iteration",
    "elapsed_s",
    "energy_sigma_to_0_eV",
    "eentro_eV",
    "stress_xx_kB",
    "pressure_kB",
    "volume_A3",
    "max_force_eV_per_A",
]


def get_parser_rows(path_file: Path, lcase: list[str]) -> dict[str, list[str]]:
    """Return whitespace-split parser rows keyed by case name."""
    dict_rows = {}
    for line in path_file.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if parts and parts[0] in lcase:
            dict_rows[parts[0]] = parts
    return dict_rows


def run_parsers(path_data: Path, lcase: list[str]) -> list[dict[str, object]]:
    """Run PostTime/PostData/PostData2 and normalize their text outputs."""
    path_data = path_data.resolve()
    lcase_present = [
        case for case in lcase if (path_data / case / "OUTCAR").is_file()
    ]

    if lcase_present:
        with NamedTemporaryFile(
                prefix="pjvasp-outcar-",
                suffix=".txt",
                delete=False) as file_temp:
            path_temp = Path(file_temp.name)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                PostTime(
                    my_path=str(path_data) + "/",
                    post_path=str(path_temp),
                ).read_OUTCAR(job_list=lcase_present)
            dict_time = get_parser_rows(path_temp, lcase_present)
            with contextlib.redirect_stdout(io.StringIO()):
                PostData(
                    my_path=str(path_data) + "/",
                    post_path=str(path_temp),
                ).read_OUTCAR(job_list=lcase_present)
            dict_data = get_parser_rows(path_temp, lcase_present)
            with contextlib.redirect_stdout(io.StringIO()):
                PostData2(
                    my_path=str(path_data) + "/",
                    post_path=str(path_temp),
                ).read_OUTCAR(job_list=lcase_present)
            dict_data2 = get_parser_rows(path_temp, lcase_present)
        finally:
            if path_temp.is_file():
                path_temp.unlink()
    else:
        dict_time = {}
        dict_data = {}
        dict_data2 = {}

    lrow = []
    for case in lcase:
        if case not in lcase_present:
            lrow.append({
                "case": case,
                "status": "missing",
                "converged": False,
                "iteration": "",
                "elapsed_s": "",
                "energy_sigma_to_0_eV": "",
                "eentro_eV": "",
                "stress_xx_kB": "",
                "pressure_kB": "",
                "volume_A3": "",
                "max_force_eV_per_A": "",
            })
            continue

        time_row = dict_time[case]
        data_row = dict_data[case]
        data2_row = dict_data2[case]
        lrow.append({
            "case": case,
            "status": "parsed",
            "converged": time_row[1] == "reachedrequiredaccuracy",
            "iteration": time_row[2].removeprefix("Iteration:"),
            "elapsed_s": float(time_row[5]),
            "energy_sigma_to_0_eV": float(data_row[1]),
            "eentro_eV": float(data_row[2]),
            "stress_xx_kB": float(data_row[5]),
            "pressure_kB": float(data2_row[2]),
            "volume_A3": float(data2_row[1]),
            "max_force_eV_per_A": float(data2_row[3]),
        })
    return lrow


def write_csv(path_output: Path, lrow: list[dict[str, object]]) -> None:
    """Write the normalized records as a deterministic CSV file."""
    path_output = path_output.resolve()
    path_output.parent.mkdir(parents=True, exist_ok=True)
    with path_output.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(lrow)


def print_table(lrow: list[dict[str, object]]) -> None:
    """Print the tutorial's compact, searchable terminal table."""
    print(
        "case   status  converged  energy_eV      pressure_kB  "
        "volume_A3  fmax_eV_A"
    )
    for row in lrow:
        if row["status"] == "missing":
            print(f"{row['case']:<6} missing False")
            continue
        print(
            f"{row['case']:<6} parsed  {str(row['converged']):<9} "
            f"{row['energy_sigma_to_0_eV']:>13.8f} "
            f"{row['pressure_kB']:>12.5f} "
            f"{row['volume_A3']:>10.2f} "
            f"{row['max_force_eV_per_A']:>10.7f}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize selected OUTCAR directories.")
    parser.add_argument(
        "path_data",
        type=Path,
        help="Directory containing case/OUTCAR subdirectories.",
    )
    parser.add_argument(
        "-cases",
        nargs="+",
        default=["0.997", "1.000"],
        help="Case directory names to parse.",
    )
    parser.add_argument(
        "-output",
        type=Path,
        default=Path("outcar_summary.csv"),
        help="CSV output path.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    lsummary = run_parsers(args.path_data, args.cases)
    assert len(lsummary) == len(args.cases)
    write_csv(args.output, lsummary)
    print_table(lsummary)
    print("wrote: " + str(args.output.resolve()))
