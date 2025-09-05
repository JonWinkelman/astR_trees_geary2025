import pandas as pd
import subprocess
from pathlib import Path
import os

def run_raxml(
    aln_fp: str | Path,
    out_dir: str | Path,
    prefix: str = "AstR",
    threads: int = 8,
    model: str = "AUTO",
    bootstraps: int = 100,
    seed: int = 12345,
    raxml_bin: str = "raxmlHPC-PTHREADS"
) -> None:
    """
    Run RAxML v8 on an alignment with ML + bootstrap.

    Parameters
    ----------
    aln_fp : str or Path
        Path to alignment file (.aln in FASTA format).
    out_dir : str or Path
        Directory to place RAxML outputs.
    prefix : str
        Name prefix for RAxML output files.
    threads : int
        Number of CPU threads to use (only with -PTHREADS binary).
    model : str
        Substitution model. "AUTO" lets RAxML pick; "LG" fixes to LG.
    bootstraps : int
        Number of bootstrap replicates.
    seed : int
        Random seed for reproducibility (used for both parsimony and rapid bootstraps).
    raxml_bin : str
        RAxML binary to call ("raxmlHPC-PTHREADS" or "raxmlHPC").
    """
    aln_fp = Path(aln_fp).resolve()
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # choose model string
    if model.upper() == "AUTO":
        model_str = "PROTGAMMAAUTO"
    elif model.upper() == "LG":
        model_str = "PROTGAMMALG"
    else:
        raise ValueError(f"Unsupported model: {model}")

    cmd = [
        raxml_bin,
        "-T", str(threads),
        "-m", model_str,
        "-p", str(seed),
        "-x", str(seed),
        "-#", str(bootstraps),
        "-f", "a",
        "-s", str(aln_fp),
        "-n", prefix,
        "-w", str(out_dir)
    ]

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"RAxML finished. Outputs in {out_dir}")


def run_iqtree_asr(alignment_file, tree_file=None, model='MFP',
                   threads=2, prefix='my_asr', output_dir=None,
                   sh_alrt_reps=1000, ufboot_reps=None, outgroups=None):
    """
    Runs IQ-TREE ancestral sequence reconstruction with SH-aLRT support.
    Optionally also runs ultrafast bootstrap (UFBoot) and supports outgroup rooting.
    
    Parameters:
    alignment_file: Path to alignment (FASTA, PHYLIP, etc.).
    tree_file: Optional Newick tree. If None, IQ-TREE infers one.
    model: Substitution model (default 'MFP').
    threads: Number of CPU threads.
    prefix: Output file prefix.
    output_dir: Directory for output files (optional).
    sh_alrt_reps: SH-aLRT replicates (default: 1000).
    ufboot_reps: UFBoot replicates (default: None → skip).
    outgroups: Outgroup taxa (comma-separated string or path to text file with one name per line).
    """
    
    # Ensure output directory exists
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        prefix_path = os.path.join(output_dir, prefix)
    else:
        prefix_path = prefix
    
    # Build base command
    command = [
        "iqtree2",
        "-s", alignment_file,
        "-m", model,
        "--ancestral",
        "-T", str(threads),
        "--alrt", str(sh_alrt_reps),
        "--prefix", prefix_path
    ]
    
    # Add bootstrap if requested
    if ufboot_reps:
        command.extend(["-B", str(ufboot_reps)])
    
    # Add tree file if provided
    if tree_file:
        command.extend(["-t", tree_file])
    
    # Handle outgroups
    if outgroups:
        if os.path.isfile(outgroups):  # If it's a file, read and join
            with open(outgroups) as f:
                outgroup_list = [line.strip() for line in f if line.strip()]
            outgroup_str = ",".join(outgroup_list)
        else:
            outgroup_str = outgroups  # Assume it's already comma-separated
        command.extend(["-o", outgroup_str])
    
    # Run IQ-TREE
    subprocess.run(command, check=True)
    print(f"IQ-TREE ASR completed. Results saved with prefix '{prefix_path}'")



import subprocess, tempfile, os
from pathlib import Path

def _detect_muscle_version(muscle_bin="muscle") -> int:
    """
    Returns major version (3 or 5). Falls back to 3 if unknown.
    """
    try:
        out = subprocess.run(
            [muscle_bin, "-version"],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False
        ).stdout.lower()
        if "muscle v5" in out or "muscle5" in out:
            return 5
        if "muscle v3" in out or "v3.8." in out:
            return 3
    except FileNotFoundError:
        raise FileNotFoundError(f"Cannot find MUSCLE executable: {muscle_bin}")
    # Some builds print version only on stderr or in different format; default to 3
    return 3

def run_muscle(inp, output_fp, muscle_bin="muscle", extra_args=None):
    """
    Run MUSCLE (v3 or v5) on sequences.

    Parameters
    ----------
    inp : dict | str | Path
        If dict: {header: sequence}. If str/Path: path to FASTA file.
    output_fp : str | Path
        Path to write alignment (FASTA).
    muscle_bin : str
        Executable name/path ('muscle', 'muscle5', or full path).
    extra_args : list[str] | None
        Additional CLI args to pass through (version-specific flags are your responsibility).

    Returns
    -------
    Path to the output alignment.
    """
    output_fp = str(output_fp)
    tmp_fa = None

    # Prepare input FASTA
    if isinstance(inp, (str, Path)):
        in_fa = str(inp)
    elif isinstance(inp, dict):
        fd, tmp_fa = tempfile.mkstemp(prefix="muscle_", suffix=".fa")
        os.close(fd)
        with open(tmp_fa, "w") as fh:
            for h, seq in inp.items():
                h = str(h).strip().lstrip(">")
                seq = "".join(str(seq).split())
                fh.write(f">{h}\n{seq}\n")
        in_fa = tmp_fa
    else:
        raise TypeError("inp must be a FASTA path or a dict {header: sequence}")

    # Detect version & build command
    major = _detect_muscle_version(muscle_bin)
    if major == 5:
        cmd = [muscle_bin, "-align", in_fa, "-output", output_fp]
    else:  # v3
        cmd = [muscle_bin, "-in", in_fa, "-out", output_fp]

    if extra_args:
        cmd.extend(extra_args)

    # Run
    try:
        print("Running:", " ".join(cmd))
        res = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if res.stdout:
            print(res.stdout.strip())
    except subprocess.CalledProcessError as e:
        msg = e.stdout or ""
        raise RuntimeError(f"MUSCLE failed (exit {e.returncode}). Output:\n{msg}") from e
    finally:
        if tmp_fa and os.path.exists(tmp_fa):
            os.remove(tmp_fa)

    return output_fp

