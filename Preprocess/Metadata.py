#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def metadata_check(DATA_DIR):

    INPUT = Path(DATA_DIR) / "Metadata_for_M2_Biomede.txt"
    COUNTS = Path(DATA_DIR) / "merged_samples_SYMBOL_raw_count.tsv"
    OUT = Path("metadata_qc")
    FIG = OUT / "figures"

    for folder in ["clinical", "experimental", "molecular", "groups"]:
        (FIG / folder).mkdir(parents=True, exist_ok=True)

    ID = ["files", "sample", "num"]
    GROUP = "DMG_subgroup"
    NUMERIC = ["age", "vaf_h3_K27M"]

    CLINICAL = [
        "Sexe", "Country", "survival_cat",
        "Metastases_at_diagnosis",
        "Metastases_at_progression",
        "Localisation", "treat", "death"
    ]

    EXPERIMENTAL = [
        "data_omics", "kit",
        "ihc_loss_H3_K27me3", "ihc_mut_K27M_H33"
    ]

    STATUS = ["H3", "TP53_status", "ACVR1_status", "PPM1D_status"]

    metadata = pd.read_csv(INPUT, sep="\t", dtype=str)
    metadata.columns = metadata.columns.str.strip()

    if metadata.empty or metadata.columns.duplicated().any():
        raise ValueError("Empty metadata or duplicated column names.")

    # 1. GENERAL OVERVIEW

    overview = pd.DataFrame({
        "column": metadata.columns,
        "missing_count": metadata.isna().sum().values,
        "missing_percentage": (metadata.isna().mean().values * 100).round(2),
        "unique_values": metadata.nunique().values,
        "constant": metadata.nunique().values <= 1
    })

    overview.to_csv(OUT / "metadata_overview.tsv", sep="\t", index=False)

    # 2. MISSING VALUES

    missing = overview[[
        "column", "missing_count", "missing_percentage"
    ]].sort_values("missing_percentage", ascending=False)

    missing.to_csv(OUT / "missing_values.tsv", sep="\t", index=False)

    # 3. VARIABLE CLASSIFICATION

    numeric_cols = [c for c in NUMERIC if c in metadata.columns]

    for col in numeric_cols:
        metadata[col] = pd.to_numeric(metadata[col], errors="coerce")

    categorical_cols = [
        c for c in metadata.columns
        if c not in ID + numeric_cols
    ]

    categorical_summary = []
    categorical_frequencies = []

    for col in categorical_cols:
        values = metadata[col].dropna()
        frequencies = values.value_counts()

        categorical_summary.append({
            "variable": col,
            "categories": values.nunique(),
            "missing": metadata[col].isna().sum()
        })

        for category, count in frequencies.items():
            categorical_frequencies.append({
                "variable": col,
                "category": category,
                "count": count,
                "percentage": count / len(values) * 100
            })

    pd.DataFrame(categorical_summary).to_csv(
        OUT / "categorical_variables.tsv", sep="\t", index=False
    )

    pd.DataFrame(categorical_frequencies).to_csv(
        OUT / "categorical_frequencies.tsv", sep="\t", index=False
    )

    # 4. NUMERICAL STATISTICS

    numerical_stats = []

    for col in numeric_cols:
        x = metadata[col]

        numerical_stats.append({
            "variable": col,
            "n": x.count(),
            "mean": x.mean(),
            "median": x.median(),
            "std": x.std(),
            "min": x.min(),
            "q1": x.quantile(0.25),
            "q3": x.quantile(0.75),
            "max": x.max()
        })

    pd.DataFrame(numerical_stats).to_csv(
        OUT / "numerical_variables.tsv", sep="\t", index=False
    )

    # 5. CONSTANT VARIABLES

    constants = metadata.columns[
        metadata.nunique(dropna=True) <= 1
    ].tolist()

    pd.DataFrame({"constant_variable": constants}).to_csv(
        OUT / "constant_variables.tsv", sep="\t", index=False
    )

    # 6. DUPLICATED IDENTIFIERS

    duplicates = []

    for col in ID:
        if col in metadata.columns:
            duplicates.append({
                "identifier": col,
                "duplicated_rows": metadata[col].duplicated(keep=False).sum(),
                "unique_values": metadata[col].nunique(),
                "missing": metadata[col].isna().sum()
            })

    pd.DataFrame(duplicates).to_csv(
        OUT / "duplicated_identifiers.tsv", sep="\t", index=False
    )

    # 7. SAMPLE CONSISTENCY

    expression_samples = pd.read_csv(
        COUNTS, sep="\t", nrows=0
    ).columns[1:].tolist()

    consistency = []

    if "sample" in metadata.columns:
        sample_ids = set(metadata["sample"].dropna())
        expression_ids = set(expression_samples)

        consistency.append({
            "metadata_samples": len(sample_ids),
            "expression_samples": len(expression_ids),
            "missing_in_expression": len(sample_ids - expression_ids),
            "missing_in_metadata": len(expression_ids - sample_ids),
            "match": sample_ids == expression_ids
        })

    pd.DataFrame(consistency).to_csv(
        OUT / "sample_consistency.tsv", sep="\t", index=False
    )

    # 8. GROUP COMPOSITION

    if GROUP in metadata.columns:
        group_counts = metadata[GROUP].value_counts(dropna=False)

        group_counts.rename_axis("group").reset_index(
            name="sample_count"
        ).to_csv(OUT / "group_composition.tsv", sep="\t", index=False)

    # 9. MUTATION ANALYSIS

    excluded = ID + numeric_cols + [GROUP] + STATUS + CLINICAL + EXPERIMENTAL

    mutation_cols = [
        c for c in metadata.columns
        if c not in excluded
        and set(metadata[c].dropna().str.upper().unique()) <= {"MUT", "WT"}
        and metadata[c].notna().any()
    ]

    for col in mutation_cols:
        metadata[col] = metadata[col].str.strip().str.upper()

    mutation_stats = []

    for col in mutation_cols:
        x = metadata[col]
        observed = x.isin(["MUT", "WT"]).sum()
        mut = (x == "MUT").sum()
        wt = (x == "WT").sum()

        mutation_stats.append({
            "gene": col,
            "MUT": mut,
            "WT": wt,
            "NA": x.isna().sum(),
            "observed": observed,
            "MUT_percentage": mut / observed * 100 if observed else np.nan
        })

    mutation_df = pd.DataFrame(mutation_stats)

    mutation_df.to_csv(
        OUT / "mutation_summary.tsv", sep="\t", index=False
    )

    if GROUP in metadata.columns:
        mutation_group = []

        for group, subset in metadata.groupby(GROUP, dropna=True):
            for col in mutation_cols:
                x = subset[col]
                observed = x.isin(["MUT", "WT"]).sum()
                mut = (x == "MUT").sum()

                mutation_group.append({
                    "group": group,
                    "gene": col,
                    "MUT": mut,
                    "WT": (x == "WT").sum(),
                    "NA": x.isna().sum(),
                    "MUT_percentage": mut / observed * 100 if observed else np.nan
                })

        pd.DataFrame(mutation_group).to_csv(
            OUT / "mutation_by_group.tsv", sep="\t", index=False
        )

    # 10. GENERIC FIGURE FUNCTIONS

    def bar_plot(data, columns, folder):
        for col in columns:
            if col not in data.columns:
                continue

            values = data[col].fillna("Missing").astype(str).value_counts()

            if values.empty:
                continue

            ax = values.plot.barh(figsize=(7, max(3, len(values) * 0.35)))
            ax.set_title(col)
            ax.set_xlabel("Number of samples")
            plt.tight_layout()
            plt.savefig(folder / f"{col}.png", dpi=300)
            plt.close()

    def numeric_plot(data, columns, folder):
        for col in columns:
            if col not in data.columns:
                continue

            x = data[col].dropna()
            if x.empty:
                continue

            plt.figure(figsize=(7, 4))
            plt.hist(x, bins="auto", edgecolor="black")
            plt.title(f"Distribution of {col}")
            plt.xlabel(col)
            plt.ylabel("Number of samples")
            plt.tight_layout()
            plt.savefig(folder / f"{col}.png", dpi=300)
            plt.close()

    def grouped_plot(data, columns, folder):
        if GROUP not in data.columns:
            return

        for col in columns:
            if col not in data.columns:
                continue

            table = pd.crosstab(
                data[GROUP], data[col].fillna("Missing"), normalize="index"
            ) * 100

            if table.empty:
                continue

            ax = table.plot.bar(stacked=True, figsize=(8, 5))
            ax.set_title(f"{col} by DMG subgroup")
            ax.set_ylabel("Percentage of samples")
            ax.legend(title=col, bbox_to_anchor=(1.02, 1))
            plt.xticks(rotation=20)
            plt.tight_layout()
            plt.savefig(folder / f"{col}_by_group.png", dpi=300)
            plt.close()

    def mutation_heatmap(data, columns, folder):
        if not columns:
            return

        matrix = data[columns].replace({
            "WT": 0, "MUT": 1
        }).fillna(2).astype(int)

        genes = matrix.columns[(matrix == 1).any()].tolist()

        if not genes:
            return

        genes = sorted(
            genes,
            key=lambda g: (matrix[g] == 1).sum(),
            reverse=True
        )[:40]

        matrix = matrix[genes]

        plt.figure(figsize=(max(9, len(genes) * 0.3), 6))
        plt.imshow(matrix, aspect="auto", interpolation="nearest",
                   vmin=0, vmax=2, cmap="viridis")
        plt.xticks(range(len(genes)), genes, rotation=90, fontsize=7)
        plt.yticks(range(len(data)), [f"Sample {i+1}" for i in range(len(data))])
        plt.xlabel("Genes")
        plt.ylabel("Samples")
        plt.title("Mutation status (WT / MUT / Missing)")
        colorbar = plt.colorbar(ticks=[0, 1, 2])
        colorbar.ax.set_yticklabels(["WT", "MUT", "Missing"])
        plt.tight_layout()
        plt.savefig(folder / "mutation_heatmap.png", dpi=300)
        plt.close()

    # 11. FIGURES

    bar_plot(metadata, CLINICAL, FIG / "clinical")
    numeric_plot(metadata, ["age"], FIG / "clinical")

    bar_plot(metadata, EXPERIMENTAL, FIG / "experimental")

    bar_plot(metadata, STATUS, FIG / "molecular")
    numeric_plot(metadata, ["vaf_h3_K27M"], FIG / "molecular")
    mutation_heatmap(metadata, mutation_cols, FIG / "molecular")

    grouped_plot(metadata, CLINICAL, FIG / "groups")
    grouped_plot(metadata, EXPERIMENTAL, FIG / "groups")

    if not mutation_df.empty:
        plot_data = mutation_df[mutation_df["MUT"] > 0].sort_values(
            "MUT_percentage"
        )

        if not plot_data.empty:
            ax = plot_data.plot.barh(
                x="gene", y="MUT_percentage", legend=False,
                figsize=(8, max(4, len(plot_data) * 0.3))
            )
            ax.set_xlabel("MUT among non-missing annotations (%)")
            ax.set_ylabel("Gene")
            ax.set_title("Mutation frequency")
            plt.tight_layout()
            plt.savefig(FIG / "molecular" / "mutation_frequency.png", dpi=300)
            plt.close()

    # 12. SUMMARY REPORT

    with open(OUT / "metadata_qc_report.txt", "w") as report:
        report.write("METADATA QC REPORT\n")
        report.write(f"Samples: {metadata.shape[0]}\n")
        report.write(f"Columns: {metadata.shape[1]}\n")
        report.write(f"Missing values: {metadata.isna().sum().sum()}\n")
        report.write(f"Constant variables: {len(constants)}\n")
        report.write(f"Mutation columns: {len(mutation_cols)}\n")
        report.write(f"Numerical variables: {len(numeric_cols)}\n")
        report.write(f"Categorical variables: {len(categorical_cols)}\n")

    print("Metadata QC completed.")
    print(f"Results: {OUT.resolve()}")


if __name__ == "__main__":
    metadata_check("Data_master")
