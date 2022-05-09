# Prediction of alternative conformations using AlphaFold 2

This repository accompanies the manuscript ["Sampling alternative conformational states of transporters and receptors with AlphaFold2"](https://elifesciences.org/articles/75751) by Diego del Alamo, Davide Sala, Hassane S. Mchaourab, and Jens Meiler. The code used to generate these models can be found in `scripts/` and was derived from the closely related repository [ColabFold](https://github.com/sokrypton/ColabFold/). This repository also includes the scripts used to plot the data, which can be found in `figures/`, as well as scripts to analyze data `analyses_scripts`. Finally, a Google Colab notebook is available for use in `notebooks/`.

The model generation code does not change the underlying AlphaFold v2.0.1 prediction pipeline (multimer prediction is not currently supported - we are working on that!). Therefore, please follow the installation instructions provided by DeepMind and review the AlphaFold2 [license](https://github.com/deepmind/alphafold/blob/main/LICENSE) and [disclaimer](https://github.com/deepmind/alphafold#license-and-disclaimer) before use (additionally, please refer to the [AlphaFold FAQ](https://alphafold.ebi.ac.uk/faq) and [ColabFold FAQ](https://github.com/sokrypton/ColabFold/blob/main/README.md)). The objective of the code contained here is to provide access to otherwise hard-to-reach settings that facilitate the generation of conformationally heterogeneous models of protein structures. Genetic and/or structural databases *do not* need to be downloaded - everything is accessible through the cloud via the [MMseqs2 API](https://github.com/soedinglab/MMseqs2).

<p align="center"><img src="https://github.com/delalamo/af2_conformations/blob/a1642c8ae1dd2e7af2c3efd06bafc569512655d0/figures/header/fig1_header.png" height="250"/></p>

*De novo* prediction of protein structures in multiple alternative conformations can usually be achieved for proteins that are absent from the AlphaFold2 training set, i.e. their structures were not determined prior to 30 April 2018. Predicting multiple conformations of proteins in the training set is, in our experience, sometimes but not usually possible.

We recommend sampling across several MSA depths. When MSAs are too shallow, the proteins are totally misfolded, whereas when they are too deep the models are conformationally uniform. The "Goldilocks range" of MSA depths that achieve the maximum number of correctly folded, but structurally diverse models seems to differ from protein to protein; in our experience, they appear to correlate with the number of amino acids. In any case, initial guesses for MSA depths can range 32-128 sequences for proteins absent from the training set and 8-64 sequences for proteins in the training set (note that this is much less than the 1000-5000 sequences that are used by AlphaFold2 by default). Once generated, these models can be analyzed using any dimensionality reduction and/or clustering algorithm; in our study we use PCA and focus mainly on the models at either extreme.

### How to use the code in this repository

Before importing the code contained in the `scripts/` folder, the user needs to install the AlphaFold source code and download the parameters to a directory named `params/`. Additional Python modules that must be installed include [Numpy](https://numpy.org/), [Requests](https://docs.python-requests.org/en/latest/), and [Logging](https://abseil.io/docs/python/guides/logging).

The scripts can be imported and used out-of-the-box to fetch multiple sequence alignments and/or templates of interest:

```python
from af2_conformations.scripts import mmseqs2

# Jobname for reference
jobname = 'T4_lysozyme'

# Amino acid sequence. Whitespace and inappropriate characters are automatically removed
sequence = ("MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNCNGVIT"
            "KDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRCALINMVFQMGETGVAGFTNSL"
            "RMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL" )
            
# PDB IDs, written uppercase with chain ID specified
pdbs = ["6LB8_A",
        "6LB8_C",
        "4PK0_A",
        "6FW2_A"]

# Initializes the Runner object that queries the MMSeqs2 server
mmseqs2_runner = mmseqs2.MMSeqs2Runner( jobname, sequence )

# Fetches the data and saves to the appropriate directory
a3m_lines, template_path = mmseqs2_runner.run_job( templates = pdbs )
```

The following code then runs a prediction without templates. Note that the `max_msa_clusters` and `max_extra_msa` options can be provided to reduce the size of the multiple sequence alignment. If these are not provided, the networks default values will be used. Additional options allow the number of recycles, as well as the number of loops through the recurrent Structure Module, to be specified.

```python
from af2_conformations.scripts import predict

predict.predict_structure_no_templates( sequence, "out.pdb",
         a3m_lines, model_id = 1, max_msa_clusters = 16,
         max_extra_msa = 32, max_recycles = 1, n_struct_module_repeats = 8 )
```

To run a prediction with templates:

```python
predict.predict_structure_from_templates( sequence, "out.pdb",
        a3m_lines, template_path = template_path,
        model_id = 1, max_msa_clusters = 16, max_extra_msa = 32,
        max_recycles = 1, n_struct_module_repeats = 8 )
```

There is also functionality to introduce mutations (e.g. alanines) across the entire MSA to remove the evolutionary evidence for specific interactions (see [here](https://www.biorxiv.org/content/10.1101/2021.11.29.470469v1) and [here](https://twitter.com/sokrypton/status/1464748132852547591) on why you would want to do this). This can be achieved as follows:

```python
# Define the mutations and introduce into the sequence and MSA
residues = [ 41,42,45,46,56,59,60,63,281,282,285,286,403,407 ]
muts = { r: "A" for r in residues }
mutated_msa = util.mutate_msa( a3m_lines, muts )
```

### Known issues

Here is a shortlist of known problems that we are currently working on:
* The MMSeqs2 server queries the PDB70, rather than the full PDB. This can cause some structures to be missed if their sequences are nearly identical to those of other PDB files.
* Multimer prediction is not currently supported.
* Custom MSAs are not currently supported.

If you find any other issues please let us know in the "issues" tab above.

### Citation

If the code in this repository has helped your scientific project, please consider citing our preprint:

```bibtex
@article {10.7554/eLife.75751,
article_type = {journal},
title = {Sampling alternative conformational states of transporters and receptors with AlphaFold2},
author = {del Alamo, Diego and Sala, Davide and Mchaourab, Hassane S and Meiler, Jens},
editor = {Robertson, Janice L and Swartz, Kenton J and Robertson, Janice L},
volume = 11,
year = 2022,
month = {mar},
pub_date = {2022-03-03},
pages = {e75751},
citation = {eLife 2022;11:e75751},
doi = {10.7554/eLife.75751},
url = {https://doi.org/10.7554/eLife.75751},
journal = {eLife},
issn = {2050-084X},
publisher = {eLife Sciences Publications, Ltd},
}
```
