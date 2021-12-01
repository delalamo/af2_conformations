import os
import numpy as np

from typing import List, NoReturn

from alphafold.data import pipeline
from alphafold.data import templates
from alphafold.data.tools import hhsearch

def mk_mock_template(
    seq: str
  ) -> dict:
  
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
  lseq = len( seq )

  # Since alphafold's model requires a template input
  # We create a blank example w/ zero input, confidence -1
  aatypes = np.array(
      templates.residue_constants.sequence_to_onehot(
          "-" * lseq,
          templates.residue_constants.HHBLITS_AA_TO_ID
        )
    )

  return {
      'template_all_atom_positions': np.zeros( ( lseq, lentype, 3 ) )[ None ],
      'template_all_atom_masks': np.zeros( ( lseq, lentype ) )[ None ],
      'template_sequence': [ f'none'.encode() ],
      'template_aatype': aatypes[ None ],
      'template_confidence_scores': np.full( lseq, -1 )[ None ],
      'template_domain_names': [ f'none'.encode() ],
      'template_release_date': [ f'none'.encode() ]
  }

###############################

def mk_template(
    seq: str,
    a3m_lines = str,
    path = str
  ) -> dict:
  
  r""" Parses templates into features

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
      binary_path="hhsearch",
      databases=[ f"{ path }/pdb70" ]
  ).query( a3m_lines )

  return templates.HhsearchHitFeaturizer(
      mmcif_dir=path,
      max_template_date="2100-01-01",
      max_hits=20,
      kalign_binary_path="kalign",
      release_dates_path=None,
      obsolete_pdbs_path=None
    ).get_templates(
      query_sequence=seq,
      hits=pipeline.parsers.parse_hhr( result )
    )

###############################

def setup_features(
    seq: str,
    a3m_lines: list,
    tfeatures_in: dict
  ) -> dict:
  
  r""" Set up features for alphafold

  Parameters
  ----------
  seq : Sequence (string)
  a3m_lines : Sequence alignment lines
  tfeatures_in : Template features

  Returns
  ----------
  Alphafold features object

  """

  msa = pipeline.parsers.parse_a3m( a3m_lines )
  return {
      **pipeline.make_sequence_features(
            sequence = seq,
            description = "none",
            num_res = len( seq )
          ), **pipeline.make_msa_features(
            msas = [ msa ]
          ), **tfeatures_in
        }
