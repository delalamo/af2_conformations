# Prediction of alternative conformations using AlphaFold 2

This repository accompanies the manuscript "Shallow multiple sequence alignments allow AlphaFold2 to predict alternative conformations of transporters and GPCRs" by Diego del Alamo, Davide Sala, Hassane S. Mchaourab, and Jens Meiler. The code used to generate these models can be found in `af2_scripts/` and was derived from the closely related repository [ColabFold](https://github.com/sokrypton/ColabFold/). This repository also includes the scripts used to plot the data, which can be found in `figures/`. Finally, a Google Colab notebook is available for use in `notebooks/`.

The model generation code does not change the underlying AlphaFold v2.0.1 prediction pipeline. Therefore, please follow the installation instructions provided by DeepMind and review the AlphaFold2 [license](https://github.com/deepmind/alphafold/blob/main/LICENSE) and [disclaimer](https://github.com/deepmind/alphafold#license-and-disclaimer) before use (additionally, please refer to the [AlphaFold FAQ](https://alphafold.ebi.ac.uk/faq) and [ColabFold FAQ](https://github.com/sokrypton/ColabFold/blob/main/README.md)). The objective of the code contained here is to provide access to otherwise hard-to-reach settings that facilitate the generation of conformationally heterogeneous models of protein structures. Genetic and/or structural databases *do not* need to be downloaded - everything is accessible through the cloud via the [MMseqs2 API](https://github.com/soedinglab/MMseqs2).

<show figure>

## Predicting alternative conformations of protein structures

*De novo* prediction of protein structures in multiple alternative conformations can usually be achieved for proteins that are absent from the AlphaFold2 training set, i.e. their structures were not determined prior to 30 April 2018. Predicting multiple conformations of proteins in the training set is, in our experience, sometimes but not usually possible.

We recommend sampling across several MSA depths. When MSAs are too shallow, the proteins are totally misfolded, whereas when they are too deep the models are conformationally uniform. The "Goldilocks range" of MSA depths that achieve the maximum nmber of correctly folded,but structurally diverse models seems to differ from protein to protein; in our experience, they appear to correlate with the number of amino acids. In any case, initial guesses for MSA depths can range 32-128 sequences for proteins absent from the training set and 8-64 sequences for proteins in the training set (note that this is much less than the 1000-5000 sequences that are used by AlphaFold2 by default). Once generated, these models can be analyzed using any dimensionality reduction and/or clustering algorithm; in our study we use PCA and focus mainly on the endpoints.