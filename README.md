# README — Build an Acinetobacter Gene Tree and Relabel for iTOL

This script downloads selected **NCBI** genomes, runs **OrthoFinder** to find orthology/HOGs, extracts **AstR** ortholog protein sequences, aligns them with **MUSCLE**, infers a **RAxML** tree, then **relabels leaf names** for clean visualization in **iTOL**.  
It uses helper utilities from `jw_utils` and `orthofinder_utils`.  

Biopython’s `Phylo` handles Newick I/O; iTOL accepts Newick uploads and annotation datasets. ([Biopython Tutorial](https://biopython.org/docs/latest/Tutorial/chapter_phylo.html))

---

## What It Does (Pipeline)

1. **Clone utilities** (`jw_utils`, `orthofinder_utils`) and import helpers.  
2. **Download NCBI metadata & files** for the assemblies in `figure_genome_assemblies`; write a summary DataFrame. Uses NCBI Datasets CLI-backed workflow and moves `GFF`/proteomes into project dirs.  
3. **Run OrthoFinder** (v2.5.5) on the proteomes to obtain orthogroups and **HOGs**.  
4. **Parse OrthoFinder outputs** with `dash_ortho_parser_GTDB` to locate the `AstR` gene’s HOG and list orthologs.  
5. **Write FASTA** for `AstR` ortholog protein sequences.  
6. **Align** with **MUSCLE** to generate a multiple sequence alignment (`*.muscle.aln`).  
7. **Infer ML tree** with **RAxML** (model = `LG`, multiple threads). Use the bipartitionsBranchLabels tree for support values.  
8. **Relabel leaves for iTOL**: parse each leaf’s accession, map via `dop_obj.accession_to_name`, write a relabeled Newick. Upload to **iTOL** for visualization.

---

## Project Layout (Key Paths Used)


