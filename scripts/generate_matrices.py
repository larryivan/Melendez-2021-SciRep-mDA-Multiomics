#!/usr/bin/env python3
"""
Build count/FPKM matrices in the same format used by this repository:
  - GSE153005_count_mat.txt
  - GSE153005_FPKM_mat.txt

Key fixes vs earlier ad-hoc script:
1) Keep only paper-relevant gene biotypes (protein_coding + selected lncRNA classes).
2) Do NOT label every non-protein gene as lncRNA.
3) Collapse pseudoautosomal duplicates (e.g., *_PAR_Y) by preferring non-PAR entries.
4) Rename SRR columns to d0_r1 ... d28_r3 for downstream DE scripts.
"""

from __future__ import annotations

import glob
import re
from pathlib import Path

import pandas as pd


# ---------- User-editable paths ----------
GTF_FILE = "ref/gencode.v29.annotation.gtf"
COUNTS_GLOB = "counts/*_counts.txt"
FPKM_GLOB_DIR = "fpkm_cufflinks/*/genes.fpkm_tracking"
FPKM_GLOB_FILE = "fpkm_cufflinks/*.genes.fpkm_tracking"

OUT_COUNT = "GSE153005_count_mat.txt"
OUT_FPKM = "GSE153005_FPKM_mat.txt"

# Map your run IDs to the paper matrix column names.
# If sample names are already d0_r1-style, they are kept as-is.
SAMPLE_RENAME = {
    "SRR12070831": "d0_r1",
    "SRR12070832": "d0_r2",
    "SRR12070833": "d0_r3",
    "SRR12070834": "d14_r1",
    "SRR12070835": "d14_r2",
    "SRR12070836": "d14_r3",
    "SRR12070837": "d14_r4",
    "SRR12070838": "d28_r1",
    "SRR12070839": "d28_r2",
    "SRR14636069": "d28_r3",
}

PREFERRED_SAMPLE_ORDER = [
    "d0_r1", "d0_r2", "d0_r3",
    "d14_r1", "d14_r2", "d14_r3", "d14_r4",
    "d28_r1", "d28_r2", "d28_r3",
]

# Matches the matrix in data/GSE153005_*_mat.txt
KEEP_SPECIFIC_TYPES = {
    "protein_coding",
    "lincRNA",
    "antisense",
    "TEC",
    "sense_intronic",
    "processed_transcript",
    "sense_overlapping",
    "non_coding",
    "macro_lncRNA",
}


def canonical_gene_id(gene_id: str) -> str:
    # Remove pseudoautosomal suffix so ENSG..._PAR_Y collapses to ENSG...
    return gene_id.replace("_PAR_Y", "")


def parse_gtf_attributes(attr_text: str) -> dict:
    out = {}
    for chunk in attr_text.strip().split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        m = re.match(r'^(\S+)\s+"([^"]+)"$', chunk)
        if m:
            out[m.group(1)] = m.group(2)
    return out


def broad_type_from_specific(specific_type: str) -> str | None:
    if specific_type == "protein_coding":
        return "protein_coding"
    if specific_type in KEEP_SPECIFIC_TYPES:
        return "lncRNA"
    return None


def normalize_sample_name(sample_name: str) -> str:
    sample_name = sample_name.strip()
    if sample_name in SAMPLE_RENAME:
        return SAMPLE_RENAME[sample_name]
    return sample_name


def reorder_sample_cols(sample_cols: list[str]) -> list[str]:
    preferred = [x for x in PREFERRED_SAMPLE_ORDER if x in sample_cols]
    others = sorted([x for x in sample_cols if x not in preferred])
    return preferred + others


def collapse_pary_rows(df: pd.DataFrame, gene_col: str) -> pd.DataFrame:
    """
    Prefer non-PAR row when both ENSGxxx and ENSGxxx_PAR_Y exist.
    """
    x = df.copy()
    x["_is_par_y"] = x[gene_col].astype(str).str.contains("_PAR_Y", regex=False)
    x["geneId"] = x[gene_col].astype(str).map(canonical_gene_id)
    x = x.sort_values(["geneId", "_is_par_y"])  # non-PAR first
    x = x.drop_duplicates(subset=["geneId"], keep="first")
    x = x.drop(columns=["_is_par_y"])
    return x


def parse_gtf_annotation(gtf_file: str) -> pd.DataFrame:
    print("🚀 Parsing GENCODE GTF:", gtf_file)
    gene_anno = {}

    with open(gtf_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue

            attrs = parse_gtf_attributes(parts[8])
            raw_gene_id = attrs.get("gene_id")
            gene_name = attrs.get("gene_name")
            gene_type = attrs.get("gene_type") or attrs.get("gene_biotype")
            if not raw_gene_id or not gene_name or not gene_type:
                continue

            broad_type = broad_type_from_specific(gene_type)
            if broad_type is None:
                continue

            gene_id = canonical_gene_id(raw_gene_id)
            is_par = raw_gene_id.endswith("_PAR_Y")

            # Keep non-PAR annotation if duplicate canonical geneId appears.
            prev = gene_anno.get(gene_id)
            if prev is None or (prev["_is_par"] and not is_par):
                gene_anno[gene_id] = {
                    "geneId": gene_id,
                    "geneName": gene_name,
                    "type": broad_type,
                    "specific_type": gene_type,
                    "_is_par": is_par,
                }

    anno_df = pd.DataFrame(gene_anno.values())
    anno_df = anno_df.drop(columns=["_is_par"]).sort_values("geneId").reset_index(drop=True)
    print(f"✅ Annotation kept genes: {len(anno_df):,}")
    return anno_df


def merge_counts(anno_df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    print("🔢 Merging HTSeq count tables...")
    count_files = sorted(glob.glob(COUNTS_GLOB))
    if not count_files:
        raise FileNotFoundError(f"No count files found with pattern: {COUNTS_GLOB}")

    count_dfs = []
    for file in count_files:
        sample_raw = Path(file).name.replace("_counts.txt", "")
        sample_name = normalize_sample_name(sample_raw)

        df = pd.read_csv(file, sep="\t", header=None, names=["geneId_raw", sample_name])
        df = df[~df["geneId_raw"].astype(str).str.startswith("__")]
        df = collapse_pary_rows(df, gene_col="geneId_raw")
        df = df[["geneId", sample_name]]
        df[sample_name] = pd.to_numeric(df[sample_name], errors="coerce").fillna(0).astype(int)
        df = df.set_index("geneId")
        count_dfs.append(df)

    count_matrix = pd.concat(count_dfs, axis=1).fillna(0).astype(int).reset_index()
    sample_cols = [c for c in count_matrix.columns if c != "geneId"]
    sample_cols = reorder_sample_cols(sample_cols)

    final_count = pd.merge(anno_df, count_matrix, on="geneId", how="inner")
    final_count = final_count[["geneId", "geneName", "type", "specific_type"] + sample_cols]
    final_count = final_count.sort_values("geneId").reset_index(drop=True)

    final_count.to_csv(OUT_COUNT, sep="\t", index=False)
    print(f"✅ Saved: {OUT_COUNT} ({len(final_count):,} genes)")
    return final_count, sample_cols


def merge_fpkm(anno_df: pd.DataFrame, sample_cols_ref: list[str]) -> pd.DataFrame:
    print("📈 Merging Cufflinks FPKM tables...")

    fpkm_files = sorted(glob.glob(FPKM_GLOB_DIR))
    if not fpkm_files:
        fpkm_files = sorted(glob.glob(FPKM_GLOB_FILE))
    if not fpkm_files:
        raise FileNotFoundError(
            f"No FPKM files found with patterns: {FPKM_GLOB_DIR} or {FPKM_GLOB_FILE}"
        )

    fpkm_dfs = []
    for file in fpkm_files:
        p = Path(file)
        if p.name == "genes.fpkm_tracking":
            sample_raw = p.parent.name
        else:
            sample_raw = p.name.replace(".genes.fpkm_tracking", "")
        sample_name = normalize_sample_name(sample_raw)

        df = pd.read_csv(file, sep="\t")
        if "tracking_id" not in df.columns or "FPKM" not in df.columns:
            raise ValueError(f"{file} missing required columns: tracking_id, FPKM")

        df = df[["tracking_id", "FPKM"]].rename(columns={"tracking_id": "geneId_raw", "FPKM": sample_name})
        df = collapse_pary_rows(df, gene_col="geneId_raw")
        df = df[["geneId", sample_name]]
        df[sample_name] = pd.to_numeric(df[sample_name], errors="coerce").fillna(0.0)

        # Apply pseudocount as in provided matrix style.
        df.loc[df[sample_name] == 0, sample_name] = 0.1

        df = df.set_index("geneId")
        fpkm_dfs.append(df)

    fpkm_matrix = pd.concat(fpkm_dfs, axis=1).fillna(0.1).reset_index()
    sample_cols = [c for c in fpkm_matrix.columns if c != "geneId"]

    # Keep count matrix order when possible.
    sample_cols = [c for c in sample_cols_ref if c in sample_cols] + sorted(
        [c for c in sample_cols if c not in sample_cols_ref]
    )

    # Keep the same gene universe as annotation/count matrix even when some genes
    # are absent from cufflinks tracking output. Missing values are set to pseudocount.
    final_fpkm = pd.merge(anno_df, fpkm_matrix, on="geneId", how="left")
    final_fpkm = final_fpkm[["geneId", "geneName", "type", "specific_type"] + sample_cols]
    final_fpkm[sample_cols] = final_fpkm[sample_cols].fillna(0.1)
    final_fpkm = final_fpkm.sort_values("geneId").reset_index(drop=True)

    final_fpkm.to_csv(OUT_FPKM, sep="\t", index=False)
    print(f"✅ Saved: {OUT_FPKM} ({len(final_fpkm):,} genes)")
    return final_fpkm


def main():
    anno_df = parse_gtf_annotation(GTF_FILE)
    _, sample_cols = merge_counts(anno_df)
    merge_fpkm(anno_df, sample_cols_ref=sample_cols)
    print("🎉 Done.")


if __name__ == "__main__":
    main()
