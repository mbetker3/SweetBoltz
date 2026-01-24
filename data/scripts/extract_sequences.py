import json
import gzip
from pathlib import Path
import sys

import gemmi

# -------------------- CONFIG --------------34------

REPO_ROOT = Path(__file__).resolve().parents[2] 

RECEPTORS_MAP = REPO_ROOT / "data/sweetdata/metadata/receptors_map.json"
RECEPTORS_DIR = REPO_ROOT / "data/sweetdata/receptors"
OUTPUT_FASTA = REPO_ROOT / "data/sweetdata/work/sequences.fasta"

DEDUPLICATE = True  # set to False if you want duplicates kept

# ------------------------------------------------


def load_structure(cif_path: Path) -> gemmi.Structure:
    if cif_path.suffix == ".gz":
        with gzip.open(cif_path, "rt") as f:
            return gemmi.read_structure_string(f.read())
    else:
        return gemmi.read_structure(str(cif_path))


def extract_chain_sequence(structure: gemmi.Structure, chain_id: str) -> str:
    # optional, but valid if your gemmi has it
    if hasattr(structure, "remove_alternative_conformations"):
        structure.remove_alternative_conformations()

    for model in structure:
        for chain in model:
            if chain.name != chain_id:
                continue

            seq = []
            for res in chain.get_polymer():  # exists on Chain
                info = gemmi.find_tabulated_residue(res.name)  # ResidueInfo
                if info and info.found() and info.is_amino_acid():
                    code = info.fasta_code()  # 1-letter-ish FASTA code
                    # fasta_code() can include lower-case/space for nonstandard; skip unknowns
                    if code and code != "?":
                        seq.append(code.upper())

            if not seq:
                raise ValueError(f"Chain {chain_id} found but no amino-acid polymer residues")
            return "".join(seq)

    raise KeyError(f"Chain {chain_id} not found in structure")


def main():
    if not RECEPTORS_MAP.exists():
        sys.exit(f"Missing receptors_map.json at {RECEPTORS_MAP}")

    OUTPUT_FASTA.parent.mkdir(parents=True, exist_ok=True)

    with open(RECEPTORS_MAP) as f:
        receptors = json.load(f)

    seen = set()
    written = 0

    with open(OUTPUT_FASTA, "w") as out:
        for pdb_id, info in receptors.items():
            state = info["state"]
            cif_path = RECEPTORS_DIR / state / f"{pdb_id}-assembly1.cif.gz"

            if not cif_path.exists():
                raise FileNotFoundError(f"Missing CIF: {cif_path}")

            structure = load_structure(cif_path)

            for role in ("tas1r2_chain", "tas1r3_chain"):
                chain_id = info[role]
                seq = extract_chain_sequence(structure, chain_id)

                if DEDUPLICATE and seq in seen:
                    continue

                seen.add(seq)
                header = f">{pdb_id}_{role}_{chain_id}"
                out.write(header + "\n")
                out.write(seq + "\n")

                written += 1

    print(f"Wrote {written} sequences to {OUTPUT_FASTA}")
    if DEDUPLICATE:
        print(f"(Deduplicated to {len(seen)} unique sequences)")


if __name__ == "__main__":
    main()