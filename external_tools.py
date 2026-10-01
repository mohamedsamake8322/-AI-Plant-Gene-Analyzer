"""Optional WSL-backed wrappers for external alignment and phylogeny tools.

The application remains usable without these binaries: the internal alignment
engine is the fallback. External programs are invoked with argument lists and
without a shell so sequence content cannot become command syntax.
"""

from __future__ import annotations

import re
import shlex
import subprocess
import tempfile
from io import StringIO
from pathlib import Path

from Bio import Phylo


WSL_DISTRIBUTION = "Ubuntu"
TOOL_COMMANDS = {
    "MAFFT": "mafft",
    "MUSCLE": "muscle",
    "ClustalW": "clustalw",
    "IQ-TREE": "iqtree2",
}


def _wsl_path(path: Path) -> str:
    resolved = path.resolve()
    if not resolved.drive:
        return resolved.as_posix()
    drive = resolved.drive.rstrip(":").lower()
    portable = resolved.as_posix()[len(resolved.drive):].lstrip("/")
    return f"/mnt/{drive}/{portable}"


def _run_wsl(arguments: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["wsl.exe", "-d", WSL_DISTRIBUTION, "--"] + arguments,
        check=True,
        capture_output=True,
        text=True,
    )


def _run_wsl_shell(
    arguments: list[str],
    stdout_path: str | None = None,
) -> subprocess.CompletedProcess[str]:
    command = shlex.join(arguments)
    if stdout_path:
        command += " > " + shlex.quote(stdout_path)
    return subprocess.run(
        ["wsl.exe", "-d", WSL_DISTRIBUTION, "--", "bash", "-lc", command],
        check=True,
        capture_output=stdout_path is None,
        text=True,
    )


def tool_status() -> dict[str, bool]:
    """Return availability of the four external tools in the WSL distro."""
    status: dict[str, bool] = {}
    for label, command in TOOL_COMMANDS.items():
        try:
            _run_wsl(["bash", "-lc", f"command -v {command}"])
        except (OSError, subprocess.CalledProcessError):
            status[label] = False
        else:
            status[label] = True
    return status


def _parse_fasta(text: str) -> tuple[list[str], list[str]]:
    labels: list[str] = []
    sequences: list[str] = []
    current_label: str | None = None
    chunks: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_label is not None:
                labels.append(current_label)
                sequences.append("".join(chunks))
            current_label = line[1:].strip()
            chunks = []
        elif current_label is not None:
            chunks.append(line)
    if current_label is not None:
        labels.append(current_label)
        sequences.append("".join(chunks))
    return labels, sequences


def restore_newick_labels(newick: str, labels_by_id: dict[str, str]) -> str:
    """Replace temporary tool-safe leaf IDs with original FASTA headers."""
    tree = Phylo.read(StringIO(newick.strip()), "newick")
    unknown_ids = [
        clade.name for clade in tree.get_terminals()
        if clade.name not in labels_by_id
    ]
    if unknown_ids:
        raise RuntimeError(f"IQ-TREE returned unknown sequence IDs: {', '.join(unknown_ids)}")

    for clade in tree.get_terminals():
        clade.name = labels_by_id[clade.name]

    output = StringIO()
    Phylo.write(
        tree,
        output,
        "newick",
        format_branch_length="%1.10f",
        format_confidence="%1.0f",
    )
    serialized = output.getvalue().strip()
    return re.sub(r":0(?:\.0+)?;$", ";", serialized)


def restore_iqtree_report_labels(report: str, labels_by_id: dict[str, str]) -> str:
    """Restore temporary sequence IDs in the human-readable IQ-TREE report."""
    for sequence_id, label in labels_by_id.items():
        report = re.sub(rf"\b{re.escape(sequence_id)}\b", lambda _: label, report)
    return report


def run_external_msa(
    sequences: list[str],
    labels: list[str],
    engine: str,
) -> dict[str, object]:
    """Run MAFFT, MUSCLE, or ClustalW and return aligned FASTA data."""
    if engine not in {"MAFFT", "MUSCLE", "ClustalW"}:
        raise ValueError(f"Unsupported external MSA engine: {engine}")
    if len(sequences) < 2:
        raise ValueError("At least two sequences are required for an external MSA.")

    with tempfile.TemporaryDirectory(prefix="plant_gene_msa_") as temp_dir:
        root = Path(temp_dir)
        input_path = root / "input.fasta"
        output_path = root / "aligned.fasta"
        input_path.write_text(
            "\n".join(f">{label}\n{sequence}" for label, sequence in zip(labels, sequences)) + "\n",
            encoding="utf-8",
        )
        wsl_input = _wsl_path(input_path)
        wsl_output = _wsl_path(output_path)

        if engine == "MAFFT":
            _run_wsl_shell(["mafft", "--auto", wsl_input], stdout_path=wsl_output)
            aligned_text = output_path.read_text(encoding="utf-8")
        elif engine == "MUSCLE":
            _run_wsl_shell(["muscle", "-align", wsl_input, "-output", wsl_output])
            aligned_text = output_path.read_text(encoding="utf-8")
        else:
            _run_wsl_shell([
                "clustalw",
                f"-INFILE={wsl_input}",
                f"-OUTFILE={wsl_output}",
                "-OUTPUT=FASTA",
                "-QUIET",
            ])
            aligned_text = output_path.read_text(encoding="utf-8")

    aligned_labels, aligned_sequences = _parse_fasta(aligned_text)
    if len(aligned_sequences) != len(sequences):
        raise RuntimeError(f"{engine} returned an incomplete alignment.")
    return {
        "algorithm": engine,
        "labels": aligned_labels,
        "aligned_sequences": aligned_sequences,
    }


def run_iqtree(
    aligned_sequences: list[str],
    labels: list[str],
    bootstrap: int = 1000,
) -> dict[str, object]:
    """Run IQ-TREE ModelFinder plus ultrafast bootstrap on an alignment."""
    if len(aligned_sequences) < 3:
        raise ValueError("At least three aligned sequences are required for a tree.")
    if bootstrap < 1000:
        raise ValueError("IQ-TREE ultrafast bootstrap requires at least 1000 replicates.")

    with tempfile.TemporaryDirectory(prefix="plant_gene_iqtree_") as temp_dir:
        root = Path(temp_dir)
        input_path = root / "alignment.fasta"
        prefix_path = root / "iqtree_result"
        input_path.write_text(
            "\n".join(f">{label}\n{sequence}" for label, sequence in zip(labels, aligned_sequences)) + "\n",
            encoding="utf-8",
        )
        _run_wsl_shell([
            "iqtree2",
            "-s", _wsl_path(input_path),
            "-m", "MFP",
            "-B", str(bootstrap),
            "-T", "AUTO",
            "--prefix", _wsl_path(prefix_path),
        ])
        tree_path = Path(str(prefix_path) + ".treefile")
        report_path = Path(str(prefix_path) + ".iqtree")
        if not tree_path.exists():
            raise RuntimeError("IQ-TREE did not produce a tree file.")
        tree = tree_path.read_text(encoding="utf-8").strip()
        report = report_path.read_text(encoding="utf-8") if report_path.exists() else ""

    model_match = re.search(r"Best-fit model.*?:\s*(\S+)", report, flags=re.IGNORECASE)
    return {
        "algorithm": "IQ-TREE",
        "model": model_match.group(1) if model_match else "MFP result in report",
        "bootstrap": bootstrap,
        "newick": tree,
        "report": report,
    }