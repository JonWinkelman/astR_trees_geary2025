# README — Build an Acinetobacter Gene Tree and Relabel for iTOL

This script downloads selected **NCBI** genomes, runs **OrthoFinder** to find orthology/HOGs, extracts **AstR** ortholog protein sequences, aligns them with **MUSCLE**, infers a **RAxML** tree, then **relabels leaf names** for clean visualization in **iTOL**.  
It uses helper utilities from `jw_utils` and `orthofinder_utils`.  

---

# What It Does (Pipeline)

1. **Clone utilities** (`jw_utils`, `orthofinder_utils`) and import helpers.  
2. **Download NCBI metadata & files** for the assemblies in `figure_genome_assemblies`; write a summary DataFrame. Uses NCBI Datasets CLI-backed workflow and moves `GFF`/proteomes into project dirs.  
3. **Run OrthoFinder** (v2.5.5) on the proteomes to obtain orthogroups and **HOGs**.  
4. **Parse OrthoFinder outputs** with `dash_ortho_parser_GTDB` to locate the `AstR` gene’s HOG and list orthologs.  
5. **Write FASTA** for `AstR` ortholog protein sequences.  
6. **Align** with **MUSCLE** to generate a multiple sequence alignment (`*.muscle.aln`).  
7. **Infer ML tree** with **RAxML** (model = `LG`, multiple threads). Use the bipartitionsBranchLabels tree for support values.  
8. **Relabel leaves for iTOL**: parse each leaf’s accession, map via `dop_obj.accession_to_name`, write a relabeled Newick. Upload to **iTOL** for visualization.

---


## Requirements

- **Python 3** with **Biopython** (`Bio.Phylo`)  
- **MUSCLE** (in `PATH`) for alignment  
- **RAxML** binary compatible with the wrapper (`ef.run_raxml`)  
- **OrthoFinder** (≥ v2.4) for orthogroup/HOG inference  
- **NCBI Datasets CLI** if running full download  

> Also requires local modules:
> * `jw_utils` — itol / ncbi / FASTA helpers  
> * `orthofinder_utils` — parsers and OrthoFinder runner  
> * `external_functions` — wrappers for MUSCLE, RAxML, etc.

---

## How to Run (assumes you have conda installed)

1. Create an environment from `environment.yml`
   `conda env create -f environment.yml'
   `conda activate Geary_astR_trees` 
3. Execute the script (adjusting or removing notebook-specific lines like `get_ipython()`).  
4. Outputs you should obtain:
   - `summary_data/AssemblyAccession_to_SpeciesName.json`
   - `astR_orthologs.faa` and `*.muscle.aln`
   - `raxML_output/RAxML_bipartitionsBranchLabels.AstR`
   - `itol_annotations/RELABLE_RAxML_bipartitionsBranchLabels.astR`
5. Upload the relabeled Newick file to **iTOL** via the Upload interface. Use iTOL’s tools to style labels, branch support, etc.

---

## Notable Implementation Details

- **Biopython Phylo** used to read Newick trees; `tree.get_terminals()` gives leaf nodes.  
- **RAxML outputs**:
  * `RAxML_bestTree.*` — best topology  
  * `RAxML_bipartitions.*` — includes support values  
  * `RAxML_bipartitionsBranchLabels.*` — supports embedded as branch labels  


