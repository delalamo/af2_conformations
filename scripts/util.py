import os
import numpy as np

from typing import Dict, List, NoReturn

from alphafold.data import pipeline
from alphafold.data import templates
from alphafold.data.tools import hhsearch


def mk_mock_template(seq: str) -> dict:

    r"""Generates mock templates that will not influence prediction
    Taken from ColabFold version 62d7558c91a9809712b022faf9d91d8b183c328c

    Parameters
    ----------
    seq: Query sequence

    Returns
    ----------
    Dictionary with blank/empty/meaningless features

    """

    # Define constants
    lentype = templates.residue_constants.atom_type_num
    lseq = len(seq)

    # Since alphafold's model requires a template input
    # We create a blank example w/ zero input, confidence -1
    aatypes = np.array(
        templates.residue_constants.sequence_to_onehot(
            "-" * lseq, templates.residue_constants.HHBLITS_AA_TO_ID
        )
    )

    return {
        "template_all_atom_positions": np.zeros((lseq, lentype, 3))[None],
        "template_all_atom_masks": np.zeros((lseq, lentype))[None],
        "template_sequence": [f"none".encode()],
        "template_aatype": aatypes[None],
        "template_confidence_scores": np.full(lseq, -1)[None],
        "template_domain_names": [f"none".encode()],
        "template_release_date": [f"none".encode()],
    }


###############################


def mk_template(seq: str, a3m_lines=str, path=str) -> dict:

    r"""Parses templates into features

    Parameters
    ----------
    seq : Query sequence
    a3m_lines : Lines form MMSeqs2 alignment
    path : Path to templates fetched using MMSeqs2

    Returns
    ----------
    Dictionary with features

    """

    result = hhsearch.HHSearch(
        binary_path="hhsearch", databases=[f"{ path }/pdb70"]
    ).query(a3m_lines)

    return templates.HhsearchHitFeaturizer(
        mmcif_dir=path,
        max_template_date="2100-01-01",
        max_hits=20,
        kalign_binary_path="kalign",
        release_dates_path=None,
        obsolete_pdbs_path=None,
    ).get_templates(query_sequence=seq, hits=pipeline.parsers.parse_hhr(result))


###############################


def setup_features(seq: str, a3m_lines: list, tfeatures_in: dict) -> dict:

    r"""Set up features for alphafold

    Parameters
    ----------
    seq : Sequence (string)
    a3m_lines : Sequence alignment lines
    tfeatures_in : Template features

    Returns
    ----------
    Alphafold features object

    """

    msa = pipeline.parsers.parse_a3m(a3m_lines)
    return {
        **pipeline.make_sequence_features(
            sequence=seq, description="none", num_res=len(seq)
        ),
        **pipeline.make_msa_features(msas=[msa]),
        **tfeatures_in,
    }


def mutate_msa(
    a3m_lines: str,
    pos_res: Dict[int, str],
) -> str:
    r"""Mutates every position in an MSA to a residue of interest

    Example usage: mutate_msa( a3m_lines, { 15: "A", 155: "A" } )
    This will mutate residues 15 and 155 to alanine throughout the MSA

    Parameters
    ----------
    a3m_lines : Sequence alignment
    pos : Position to change
    target_res : Residue to mutate to

    Returns
    ----------
    Sequence alignment (as string)

    """

    for target_res in pos_res.values():
        assert len(target_res) == 1

    output = []

    # Iterate over alignment lines
    for line in a3m_lines.split("\n"):
        if line.startswith(">"):
            output.append(line)
        elif len(line) > 1:
            line = list(line)
            for pos, res in pos_res.items():
                if line[pos] in "ACDEFGHIKLMNPQRSTVWY":
                    line[pos] = res
            output.append("".join(line))
        else:
            output.append(line)
    return "\n".join(output)


def mutate(x, y):
    mutate_msa(x, y)  # Alias for brevity


def plddt_to_bfactor(filename: str, maxval: float = 100.0) -> NoReturn:
    r"""Converts a pLDDT vals to a B factor
    This equation is derived from the following publication:
    "Improved protein structure refinement guided by deep learning based
    accuracy estimation" by Hiranuma et al 2021
    https://doi.org/10.1038/s41467-021-21511-x

    Parameters
    ----------
    filename : Name of PDB file
    maxval : Set to 100 if using AF2 (or 1 if RoseTTAFold)

    Returns
    ----------
    None

    """
    pdb = Bio.PDB.PDBParser().get_structure("TEMP", filename)
    for atom in pdb.get_atoms():
        rmsf = 1.5 * np.exp(4 * (0.7 - (atom.bfactor / maxval)))
        atom.bfactor = (8.0 / 3.0) * (np.pi**2) * (rmsf**2)

    pdbio = Bio.PDB.PDBIO()
    pdbio.set_structure(pdb)
    pdbio.save(filename)
