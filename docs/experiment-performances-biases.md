# Evaluation of experiment performances and biases

MetaPepView provides a number of tools to perform validation and benchmarking of metaproteomics experiments. Two modules are dedicated to validation of metaproteomics experiments:

- The *Experimental performance* module: Visualization and benchmarking of experimental performance parameters from spectral data (*mzML*) + proteomics data.
- The *Evaluate community composition* module: Validate database quality used in DB search matching by comparing taxonomy composition from db search matching + local db annotation against taxonomy composition from peptides (DB search and de novo) annotated by [Unipept](https://unipept.ugent.be/).

![[Evaluation_modules.PNG]]
*Experiment validation modules highlighted in the sidebar*

In addition to dedicated experiment evaluation modules, MetaPepView provides a tool for investigation potential over- or underrepresentation of taxa in community compositions: the *Taxonomic drop-off*. This tool is part of the *[[taxonomy-visualization#Profiling peptide uniqueness in taxonomy group: Taxonomic drop-off|Community composition]]* module.



