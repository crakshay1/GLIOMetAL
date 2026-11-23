#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pathlib import Path

def Integrity(DATA_DIR):
    """
        This module takes a file as a parameter and checks for missing, infitinte, negative values, 
        duplicated identifiers and sample columns, sample consistency, matrix dimensions and identifier consistency.

    """

    # CONFIGURATION
    MATRICES = {
        "gene_raw_counts": DATA_DIR / "merged_samples_SYMBOL_raw_count.tsv",
        "gene_tpm": DATA_DIR / "merged_samples_SYMBOL_tpm.tsv",
        "transcript_raw_counts": DATA_DIR / "merged_samples_txID_raw_count.tsv",
        "transcript_tpm": DATA_DIR / "merged_samples_txID_tpm.tsv"
    }

    SALMON_DIR = DATA_DIR / "Data_for_M2_Biomede"
    REPORT_DIR = Path("INTEGRITY")
    REPORT_DIR.mkdir(exist_ok=True)
    results = []



    # REPORTING
    def add_result(matrix, check, status, details=""):
        results.append({
            "Matrix": matrix,
            "Check": check,
            "Status": status,
            "Details": str(details)
        })


    # MATRIX INTEGRITY CHECKS
    def check_matrix(name, filepath):
        print(f"\n{'=' * 60}")
        print(f"Checking: {name}")
        print(f"{'=' * 60}")

        if not filepath.exists():
            add_result(name, "File existence", "FAIL", "File not found")
            return None

        try:
            df = pd.read_csv(filepath, sep="\t", index_col=0)
        except Exception as e:
            add_result(name, "File reading", "FAIL", e)
            return None

        # Matrix dimensions
        n_rows, n_cols = df.shape
        print(f"Dimensions: {n_rows} features x {n_cols} samples")
        add_result(
            name, "Matrix dimensions", "PASS",
            f"{n_rows} features, {n_cols} samples"
        )

        # Missing values
        missing = int(df.isna().sum().sum())
        add_result(
            name, "Missing values",
            "PASS" if missing == 0 else "FAIL",
            f"{missing} missing values"
        )

        # Infinite values
        numeric = df.apply(pd.to_numeric, errors="coerce")
        infinite = int(np.isinf(numeric.to_numpy(dtype=float)).sum())
        add_result(
            name, "Infinite values",
            "PASS" if infinite == 0 else "FAIL",
            f"{infinite} infinite values"
        )

        # Negative values
        negative = int((numeric < 0).sum().sum())
        add_result(
            name, "Negative values",
            "PASS" if negative == 0 else "FAIL",
            f"{negative} negative values"
        )

        # Non-numeric values
        non_numeric = int(
            (numeric.isna() & df.notna()).sum().sum()
        )
        add_result(
            name, "Non-numeric values",
            "PASS" if non_numeric == 0 else "FAIL",
            f"{non_numeric} non-numeric values"
        )

        #  Duplicated feature identifiers
        duplicated_features = int(df.index.duplicated().sum())
        add_result(
            name, "Duplicated feature IDs",
            "PASS" if duplicated_features == 0 else "WARNING",
            f"{duplicated_features} duplicated identifiers"
        )

        # Duplicated sample columns
        duplicated_samples = int(df.columns.duplicated().sum())
        add_result(
            name, "Duplicated sample columns",
            "PASS" if duplicated_samples == 0 else "FAIL",
            f"{duplicated_samples} duplicated columns"
        )

        # Empty identifiers
        empty_ids = int(
            df.index.isna().sum() +
            (df.index.astype(str).str.strip() == "").sum()
        )
        add_result(
            name, "Empty feature IDs",
            "PASS" if empty_ids == 0 else "FAIL",
            f"{empty_ids} empty identifiers"
        )

        # Empty sample names
        empty_samples = int(
            df.columns.isna().sum() +
            (df.columns.astype(str).str.strip() == "").sum()
        )
        add_result(
            name, "Empty sample names",
            "PASS" if empty_samples == 0 else "FAIL",
            f"{empty_samples} empty sample names"
        )

        # Numeric conversion validation
        if non_numeric > 0:
            print(f"WARNING: {non_numeric} non-numeric entries")

        return {
            "name": name,
            "filepath": filepath,
            "data": df,
            "samples": list(df.columns),
            "features": df.index
        }


    # SALMON FILE INTEGRITY

    def check_salmon():
        quant_files = sorted(SALMON_DIR.glob("*/quant.sf"))
        print(f"\nSalmon quant.sf files: {len(quant_files)}")
        add_result(
            "Salmon",
            "Number of quant.sf files",
            "PASS" if len(quant_files) == 30 else "WARNING",
            f"{len(quant_files)} files found; expected 30"
        )

        expected_columns = {
            "Name", "Length", "EffectiveLength", "TPM", "NumReads"
        }

        for filepath in quant_files:
            sample = filepath.parent.name
            try:
                df = pd.read_csv(filepath, sep="\t")

                # Column validation
                missing_columns = expected_columns - set(df.columns)

                add_result(
                    sample,
                    "Salmon columns",
                    "PASS" if not missing_columns else "FAIL",
                    f"Missing columns: {missing_columns}"
                )

                if missing_columns:
                    continue

                # Missing values
                missing = int(df.isna().sum().sum())

                add_result(
                    sample,
                    "Salmon missing values",
                    "PASS" if missing == 0 else "FAIL",
                    f"{missing} missing values"
                )

                # Duplicated transcript IDs
                duplicated = int(df["Name"].duplicated().sum())

                add_result(
                    sample,
                    "Salmon duplicated transcript IDs",
                    "PASS" if duplicated == 0 else "WARNING",
                    f"{duplicated} duplicated IDs"
                )

                # Numeric columns
                numeric_cols = [
                    "Length", "EffectiveLength", "TPM", "NumReads"
                ]

                numeric = df[numeric_cols].apply(
                    pd.to_numeric, errors="coerce"
                )

                non_numeric = int(
                    (numeric.isna() & df[numeric_cols].notna()).sum().sum()
                )

                infinite = int(
                    np.isinf(numeric.to_numpy(dtype=float)).sum()
                )

                negative = int((numeric < 0).sum().sum())

                add_result(
                    sample,
                    "Salmon numerical integrity",
                    "PASS" if non_numeric + infinite + negative == 0
                    else "FAIL",
                    f"Non-numeric: {non_numeric}; "
                    f"Infinite: {infinite}; Negative: {negative}"
                )

            except Exception as e:
                add_result(sample, "Salmon file reading", "FAIL", e)


    # SAMPLE CONSISTENCY
    def check_sample_consistency(matrices):
        valid_matrices = [
            m for m in matrices if m is not None
        ]

        if not valid_matrices:
            return

        reference = set(valid_matrices[0]["samples"])
        for matrix in valid_matrices[1:]:
            current = set(matrix["samples"])
            missing = reference - current
            extra = current - reference
            status = "PASS" if not missing and not extra else "FAIL"

            add_result(
                matrix["name"],
                "Sample consistency",
                status,
                f"Missing: {len(missing)}; Extra: {len(extra)}"
            )

    matrices = []
    for name, filepath in MATRICES.items():
        matrix = check_matrix(name, filepath)
        matrices.append(matrix)

    check_salmon()
    check_sample_consistency(matrices)

    # Generate report
    report = pd.DataFrame(results)
    report.to_csv( REPORT_DIR / "integrity_report.csv", index=False)

    print("\n" + "=" * 60)
    print("DATA INTEGRITY SUMMARY")
    print("=" * 60)
    print(report.groupby(["Status"]).size())
    print(f"\nReport saved in: {REPORT_DIR}/integrity_report.csv")
