#!/usr/bin/env python3
import os
import sys
import shutil
import shlex
from pathlib import Path

import easygui

from Metadata import metadata_check
from Integrity import Integrity
from Quantification import quant_qc


def run_r(script, data_dir, output_dir):
    command = f"Rscript {shlex.quote(script)} {shlex.quote(str(data_dir))} {shlex.quote(str(output_dir))}"
    status = os.system(command)
    if status != 0:
        sys.exit(f"Failed: {script}")


def main():
    selected = easygui.diropenbox(title="Pls select folder por favor")
    if not selected:
        sys.exit("No folder selected.")
    data_dir = Path(selected).resolve()

    filtered_dir = Path("FILTERING").resolve()
    normalized_dir = Path("NORMALIZATION").resolve()
    descriptive_dir = Path("DESCRIPTIVE").resolve()

    metadata_check(str(data_dir))
    Integrity(data_dir)
    quant_qc(data_dir)

    run_r("Filter.r", data_dir, filtered_dir)
    run_r("Normalization.r", filtered_dir, normalized_dir)

    shutil.copy(data_dir / "merged_samples_SYMBOL_tpm.tsv", normalized_dir)
    shutil.copy(data_dir / "Metadata_for_M2_Biomede.txt", normalized_dir)

    run_r("Desc.r", normalized_dir, descriptive_dir)


if __name__ == "__main__":
    main()
