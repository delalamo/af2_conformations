from . import util
import os
import numpy as np
import random
import sys

from alphafold.common import protein
from alphafold.model import data
from alphafold.model import config
from alphafold.model import model

from typing import List, NoReturn

from absl import logging

def set_config(
    use_templates: bool,
    max_msa_clusters: int,
    max_extra_msa: int,
    max_recycles: int,
    model_id: int,
    n_struct_module_repeats: int,
    n_features_in: int,
    model_params: int = 0
  ): # -> alphafold.model.RunModel:
  
  r""" Generated Runner object for AlphaFold
  
  Parameters
  ----------
  use_templates : Whether templates are used
  max_msa_cluster : How many sequences to use in MSA
  max_extra_msa : How many extra sequences to use? Unclear
  max_recycles : Number of recycling iterations
  model_id : Which AF2 model to use
  n_struct_module_repeats : Number of passes through structure module
  n_features_in : Unclear
  model_params : Which AF2 model config to use

  Returns
  ----------
  AlphaFold RunModel object

  """  

  if model_id not in range( 1, 6 ):
    raise ValueError( "model_id must be between 1 and 5!" )

  # Match model_params to model_id
  # Sometimes we don't want to do this, for example,
  #   to reproduce output from ColabFold (which only uses models 1 and 3)

  logging.debug( "Prediction parameters:" )
  logging.debug( f"\tModel ID: { model_id }" )
  logging.debug( f"\tUsing templates: { use_templates }" )
  logging.debug( f"\tMaximum number of MSA clusters: { max_msa_clusters }" )
  logging.debug( f"\tMaximum number of extra MSA clusters: { max_extra_msa }" )
  logging.debug( f"\tMaximum number of recycling iterations: { max_recycles }" )

  name = f"model_{ model_params }_ptm"

  cfg = config.model_config( name )

  #### Provide config settings

  #### MSAs

  cfg.data.eval.num_ensemble = 1
  cfg.data.eval.max_msa_clusters = min(
      n_features_in, max_msa_clusters
    )
  cfg.data.common.max_extra_msa = max( 1, min(
      n_features_in - max_msa_clusters, max_extra_msa )
    )

  #### Recycle and number of iterations

  cfg.data.common.num_recycle = max_recycles
  cfg.model.num_recycle = max_recycles
  cfg.model.heads.structure_module.num_layer = n_struct_module_repeats

  #### Templates

  t = use_templates # for brevity

  cfg.data.common.use_templates = use_templates
  cfg.model.embeddings_and_evoformer.template.embed_torsion_angles = t
  cfg.model.embeddings_and_evoformer.template.enabled = t
  cfg.data.common.reduce_msa_clusters_by_max_templates = t
  cfg.data.eval.subsample_templates = t

  p = data.get_model_haiku_params(
      model_name=name, data_dir="."
    )

  return model.RunModel( cfg, p )

def run_one_job(
    runner: model.RunModel,
    features_in: dict,
    random_seed: int,
    outname: str
  ) -> NoReturn:
  r""" Runs one AF2 job with input parameters

  Parameters
  ----------
  runner : AlphaFold2 job runner
  features_in : Input features, including MSA and templates
  random_seed : Random seed
  outname : Name of PDB file to write

  Returns
  ----------
  None

  """

  # Do one last bit of processing
  features = runner.process_features(
      features_in,
      random_seed=random_seed
    )

  # Generate the model
  result = runner.predict( features, random_seed=random_seed )
  pred = protein.from_prediction( features, result )

  # Write to file
  to_pdb(
      outname,
      pred,
      result[ 'plddt' ],
      features_in[ 'residue_index' ]
    )

def predict_structure_from_templates(
    seq: str,
    outname: str,
    a3m_lines: str,
    templates: List[ str ],
    template_path: str,
    model_id: int = -1,
    model_params: int = -1,
    random_seed: int = -1,
    max_msa_clusters: int = 512,
    max_extra_msa: int = 1024,
    max_recycles: int = 3,
    n_struct_module_repeats: int = 8
  ) -> NoReturn:
  
  r""" Predicts the structure.

  Parameters
  ----------
  seq : Sequence
  outname : Name of output PDB
  a3m_lines : String of entire alignment
  templates : Choice of templates to use (format: WXYZ_A)
  template_paths : Where to locate templates
  model_id : Which AF2 model to run (must be 1 or 2 for templates)
  model_params : Which parameters to provide to AF2 model
  random_seed : Random seed
  max_msa_clusters : Number of sequences to use
  max_extra_msa : ???
  max_recycles : Number of iterations through AF2
  n_struct_module_repeats : Number of passes through structural refinement
  move_prefix : Prefix for temporary files (deleted after fxn completion)

  Returns
  ----------
  None

  """

  if random_seed == -1:
    random_seed = random.randrange( sys.maxsize )

  if model_id not in ( 1, 2 ):
    model_id = random.randint( 1, 2 )

  if model_params not in ( 1, 2 ):
    model_params = random.randint( 1, 2 )

  # Assemble the dictionary of input features
  features_in = util.setup_features(
      seq,
      a3m_lines,
      util.mk_template(
          seq,
          a3m_lines,
          template_path
        ).features
    )

  # Run the models
  model_runner = set_config(
      True,
      max_msa_clusters,
      max_extra_msa,
      max_recycles,
      model_id,
      n_struct_module_repeats,
      len( features_in[ "msa" ] ),
      model_params=model_params
    )

  run_one_job(
      model_runner,
      features_in,
      random_seed,
      outname
    )

  del model_runner

def predict_structure_no_templates(
    seq: str,
    outname: str,
    a3m_lines: str,
    model_id: int = -1,
    model_params: int = -1,
    random_seed: int = -1,
    max_msa_clusters: int = 512,
    max_extra_msa: int = 1024,
    max_recycles: int = 3,
    n_struct_module_repeats: int = 8
  ) -> NoReturn:
  
  r""" Predicts the structure.

  Parameters
  ----------
  seq : Sequence
  outname : Name of output PDB
  a3m_lines : String of entire alignment
  model_id : Which AF2 model to run (must be 1 or 2 for templates)
  random_seed : Random seed
  max_msa_clusters : Number of sequences to use
  max_extra_msa : ???
  max_recycles : Number of iterations through AF2
  n_struct_module_repeats : Number of passes through structural refinement
  
  Returns
  ----------
  None

  """

  # Set AF2 model details
  if model_id not in range( 1, 6 ):
    model_id = random.randint( 1, 5 )

  if model_params not in range( 1, 6 ):
    model_params = model_id

  if random_seed == -1:
    random_seed = random.randrange( sys.maxsize )

  features_in = util.setup_features(
      seq,
      a3m_lines,
      util.mk_mock_template( seq )
    )

  model_runner = set_config(
      False,
      max_msa_clusters,
      max_extra_msa,
      max_recycles,
      model_id,
      n_struct_module_repeats,
      len( features_in[ "msa" ] ),
      model_params=model_params
    )

  run_one_job(
      model_runner,
      features_in,
      random_seed,
      outname
    )

  del model_runner

def to_pdb(
    outname,
    pred, # type unknown but check?
    plddts, # type unknown but check?
    res_idx
  ) -> NoReturn:
  
  r""" Writes unrelaxed PDB to file

  Parameters
  ----------
  outname : Name of output PDB
  pred : Prediction to write to PDB
  plddts : Predicted errors
  res_idx : Residues to print (default=all)
  
  Returns
  ----------
  None

  """
 
  with open( outname, 'w' ) as outfile:
    outfile.write( protein.to_pdb( pred ) )

  with open( f"b_{ outname }", "w" ) as outfile:
    for line in open( outname, "r" ).readlines():
      if line[ 0:6 ] == "ATOM  ":
        seq_id = int( line[ 22:26 ].strip() ) - 1
        seq_id = np.where( res_idx == seq_id )[ 0 ][ 0 ]
        outfile.write( "{}A{}{:6.2f}{}".format(
            line[ :21 ],
            line[ 22:60 ],
            plddts[ seq_id ],
            line[ 66: ] )
        )

  os.rename( f"b_{ outname }", outname )
