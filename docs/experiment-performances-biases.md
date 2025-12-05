# Evaluation of experiment performances and biases

Meta-PepView provides a number of tools to perform evaluation and benchmarking of metaproteomics experiments. Two modules are dedicated to evaluation of metaproteomics experiments:

- The *Experiment performance evaluation* module: Visualization and benchmarking of experimental performance parameters from spectral data (*mzML*) + proteomics data. This module processes and visualizes the *performance evaluation dataset* stored in the project.
- The *Community composition evaluation* module: Validate database quality used in DB search matching by comparing taxonomy composition from DB search matching + local db annotation against taxonomy composition from peptides (DB search and *de novo*) annotated by [Unipept](https://unipept.ugent.be/). Similarly to the *Community composition* module, this module visualizes the *compositions and functions dataset* stored in the project.

![[Evaluation_modules.PNG]]
*Experiment evaluation modules highlighted in the sidebar*

In addition to dedicated experiment evaluation modules, meta-PepView provides a tool for investigation potential over- or underrepresentation of taxa in community compositions: the *Taxonomic drop-off*. This tool is part of the *[[taxonomy-visualization#Taxonomic drop-off: Profiling peptide uniqueness in taxonomy group|Community composition]]* module.



