import re
import numpy as np
import pandas as pd
from pathlib import Path
from jw_utils import parse_fasta as pfa
from jw_utils import parse_gff as pgf



def make_merged_df(gff_fp, proteome_fp):
    """Merge two dfs based on exact matches of start/end in fasta and gff file

    for mapping fasta IDs to gff IDs by start/stop when fasta/gff IDs don't match.
    """

    gff_df = pgf.make_simple_annot_df(gff_fp, start_end=True, contig=True).reset_index()
    fasta_d = pfa.get_seq_dict(str(proteome_fp), use_full_header=True)
    fasta_df = pd.DataFrame(list(fasta_d.keys()))[0].str.split('#').apply(pd.Series)
    fasta_df.columns = ['ID', 'start', 'end', 'strand', '']
    fasta_df['strand'] = fasta_df['strand'].apply(lambda x: '+' if x.strip() == '1' else '-')
    fasta_df['contig']=fasta_df['ID'].str.split('_')[:-1].apply(lambda x: '_'.join(x[:-1]))
    return match_intervals(fasta_df,gff_df)




def match_intervals(
    fasta_df: pd.DataFrame,
    gff_df: pd.DataFrame,
    keep_all: bool = False,
    fasta_suffix: str = "_fasta",
    gff_suffix: str = "_gff",
) -> pd.DataFrame:
    """
    Match intervals from fasta_df to gff_df.

    Assumptions:
      - Both dataframes have columns 'start' and 'end' (1-based inclusive).
      - If present, columns 'contig' and/or 'strand' will be used to constrain matches.

    Parameters
    ----------
    fasta_df : pd.DataFrame
    gff_df   : pd.DataFrame
   
    fasta_suffix, gff_suffix : str
        Suffixes for overlapping column names in the merged output.

    Returns
    -------
    pd.DataFrame
        Rows from fasta_df with matched gff_df columns appended.
        If no match is found for a FASTA row, it will be absent
    """
    # Work on copies; ensure integer coords
    f = fasta_df.copy()
    g = gff_df.copy()
    for c in ("start", "end"):
        if c in f: f[c] = f[c].astype(int)
        if c in g: g[c] = g[c].astype(int)

    # Decide grouping keys that must match
    keys = []
    for k in ("contig", "strand"):
        if k in f.columns and k in g.columns:
            keys.append(k)

 
    # exact start/end (and keys) match
    return pd.merge(
        f, g, how="inner",
        left_on=keys + ["start", "end"],
        right_on=keys + ["start", "end"],
        suffixes=(fasta_suffix, gff_suffix)
    )



# regex to detect protein-like IDs (avoid putting these into gene=)
PROT_LIKE = re.compile(r"^(WP_\d+(?:\.\d+)?|[A-Z]{3,}\d+(?:\.\d+)?|[A-Z]{1,6}\d{3,})$")

def write_gff_from_df_fast(df: pd.DataFrame, out_path: str | Path, source: str = "GTDBmap") -> None:
    """
    Vectorized GFF3 writer (no iterrows).
    
    Parameters
    ----------
    df : pandas.DataFrame
        Must contain:
          - contig, start, end, strand
          - gene_ID
          - ID (protein identifier, e.g. WP_ accession)
          - product
        Optional:
          - common_name
          - locus_tag
    out_path : str or Path
        Path to write GFF3 file.
    source : str
        Value for the 'source' column (default: 'GTDBmap').
    
    Behavior
    --------
    - Writes one 'gene' feature per (contig, strand, gene_ID).
    - Writes one 'CDS' feature per row.
    - Attributes are concatenated with ';' and missing values omitted.
    """

    out_path = Path(out_path)
    df = df.copy()

    # normalize coords and strand
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    if "strand" in df.columns:
        df["strand"] = (df["strand"]
                        .replace({"-1": "-", -1: "-", "1": "+", 1: "+"})
                        .fillna("+")
                        .astype(str))
    else:
        df["strand"] = "+"

    if "locus_tag" not in df.columns:
        df["locus_tag"] = "."

    if "common_name" not in df.columns:
        df["common_name"] = "."

    # ---------- GENE table ----------
    spans = (df.groupby(["contig", "strand", "gene_ID"], as_index=False)
               .agg(gstart=("start", "min"), gend=("end", "max")))

    reps = (df.sort_values(["contig", "strand", "start"])
              .drop_duplicates(["contig", "strand", "gene_ID"])[
                  ["contig", "strand", "gene_ID", "locus_tag", "common_name"]
              ])

    genes = spans.merge(reps, on=["contig", "strand", "gene_ID"], how="left")

    # blank out protein-like common_names
    mask_prot_like = genes["common_name"].astype(str).str.match(PROT_LIKE, na=False)
    genes.loc[mask_prot_like, "common_name"] = "."

    g_attr_parts = pd.DataFrame({
        "ID":        "ID=" + genes["gene_ID"].astype(str),
        "locus_tag": np.where(genes["locus_tag"].ne("."), "locus_tag=" + genes["locus_tag"].astype(str), pd.NA),
        "Name":      np.where(genes["common_name"].ne("."), "Name=" + genes["common_name"].astype(str), pd.NA),
        "gene":      np.where(genes["common_name"].ne("."), "gene=" + genes["common_name"].astype(str), pd.NA)
    })

    g_attr = (g_attr_parts.stack(future_stack=True)
                               .dropna()
                               .astype(str)
                               .groupby(level=0, sort=False)
                               .agg(";".join))

    gene_tbl = pd.DataFrame({
        "seqid":  genes["contig"].astype(str),
        "source": source,
        "type":   "gene",
        "start":  genes["gstart"].astype(int),
        "end":    genes["gend"].astype(int),
        "score":  ".",
        "strand": genes["strand"].astype(str),
        "phase":  ".",
        "attributes": g_attr.values
    })

    # ---------- CDS table ----------
    cds_gene_sym = df["common_name"].astype(str)
    cds_gene_sym = cds_gene_sym.mask(cds_gene_sym.str.match(PROT_LIKE, na=False), ".")

    coords_id = (df["contig"].astype(str) + "_" +
                 df["start"].astype(str) + "_" +
                 df["end"].astype(str) + "_" +
                 df["strand"].astype(str))

    cds_id = df["ID"].astype(str).where(df["ID"].notna() & df["ID"].ne("."), coords_id)

    c_attr_parts = pd.DataFrame({
        "ID":        "ID=" + cds_id,
        "Parent":    "Parent=" + df["gene_ID"].astype(str),
        "locus_tag": np.where(df["locus_tag"].ne("."), "locus_tag=" + df["locus_tag"].astype(str), pd.NA),
        "protein":   np.where(df["ID"].ne("."), "protein_id=" + df["ID"].astype(str), pd.NA),
        "Name":      np.where(df["ID"].ne("."), "Name=" + df["ID"].astype(str), pd.NA),
        "gene":      np.where(cds_gene_sym.ne("."), "gene=" + cds_gene_sym, pd.NA),
        "product":   np.where(df["product"].ne("."), "product=" + df["product"].astype(str), pd.NA)
    })

    c_attr = (c_attr_parts.stack(future_stack=True)
                               .dropna()
                               .astype(str)
                               .groupby(level=0, sort=False)
                               .agg(";".join))

    cds_tbl = pd.DataFrame({
        "seqid":  df["contig"].astype(str),
        "source": source,
        "type":   "CDS",
        "start":  df["start"].astype(int),
        "end":    df["end"].astype(int),
        "score":  ".",
        "strand": df["strand"].astype(str),
        "phase":  "0",
        "attributes": c_attr.values
    })

    # ---------- write ----------
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("##gff-version 3\n")
    gene_tbl.to_csv(out_path, sep="\t", header=False, index=False, mode="a")
    cds_tbl.to_csv(out_path,  sep="\t", header=False, index=False, mode="a")
