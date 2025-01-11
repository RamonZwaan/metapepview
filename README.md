
# Description
MetaPepView is a dashboard for visualization of taxonomy and protein functions from community microbiomes through metaproteomics data, as well as performance analysis of metaproteomics experiments. It is designed as a web interface that can be ran locally on desktop PC's.

Main features of the dashboard include:

- Visualize taxonomic compositions and protein function groups across microbial community samples from metaproteomics data generated from various db search or de novo search engines.
- Perform de novo taxonomy classification without a local protein-taxonomy mapping database using the Unipept API and compare compositions against db search derived taxonomy classifications.
- Analyze MS spectral performance from mzML input data and link spectral data to db search and de novo matches.
- Benchmark sample MS performance against a reference MS dataset of experiments.
# Installation / Set up
## PC requirements

- The dashboard is cross-platform and has been tested in Windows 10/11 and Linux.
- A minimum RAM of 16 GB is required for full operation of the dashboard, 32 GB is recommended.
- The dashboard can be accessed through both Firefox (Gecko browser engine) and any Chromium based browser (Chrome, Edge, Brave, etc.). **Note: the Chomium browser engine supports file import up to ~300 MB. Any larger file (common in mzML) will not be imported. Firefox does not have such a limit and multi GB files have been successfully imported.**
- Internet connection is required to download taxonomy and function datasets, and for connection to the Unipept API for de novo taxonomy classification.
## Setting up the server
MetaPepView runs on top of the dash framework in python. To access the dashboard, the server needs to be started locally on the PC. The easiest method is to start a [Docker](https://www.docker.com/) container from the DockerFile provided in the repository. For this, docker (or another OCI compliant manager) has to be installed on the system.

First, a [Docker image](https://docs.docker.com/get-started/docker-concepts/building-images/build-tag-and-publish-an-image/) needs to be built: 
> [!NOTE]
> Make sure the docker daemon is running during the command execution (start docker desktop if installed).
```Bash

$ docker build -t Metapepview . # directory of Dockerfile

```
Once the image is successfully built, start a container by running the following command:
```Bash
$ docker run -p 8050:8050 Metapepview
```
`-p` publishes port 8050 from the container to port 8050 in the host pc. This allows us to reach the dashboard in the container. We can access the dashboard by typing in `http://localhost:8050` in the web browser. Note that the dashboard is only accessible from the host PC as long as the port is closed in the firewall.

Alternatively the repository provides a requirements file to install all required python packages with [pip](https://packaging.python.org/en/latest/tutorials/installing-packages/#requirements-files). Then, the dashboard can be run directly on the host pc by executing:
```Shell
$python ./index.py
```

# Dashboard operation

## Getting started

> [!WARNING]
> The Dash framework is stateless. This means that if you refresh the web page, all data stored in memory will be wiped and a new session will be started. Make sure to export experiment data after data has been imported with *Export json*. That way, data can quickly be imported in a new session with *Import json*.

When the MetaPepView server is running and you open it in the browser at `http://localhost:8050`, you will be brought to the *Annotation* page. Here, you can import samples into a Metapepview table. This is a dataset that combines db search data, de novo peptide data and taxonomy + function classification data from multiple samples.
### Import online reference datasets.
To map peptide taxonomies against phylogenetic lineages and functions into metabolic groups or function domains, reference databases need to be imported. This is managed **on the left sidebar at the bottom** at *Database loaded*. By clicking *Fetch source*, you can download the latest version of NCBI taxonomy, GTDB taxonomy and KEGG function data from their respective sources. When using the dashboard for the first time, select the datasets for download and **create parent directories** for the data. The datasets range in the **1-5 GB** in size.

### *Annotation*: Import new sample
The dashboard provides great flexibility in importing metaproteomics dataset for for community characterization. Several formats are supported:
> [!NOTE]
> For efficient RAM handling, the dashboard supports reading of compressed datasets. For speed improvement and lower RAM usage, compress large files (>100 MB) in `.zip` or `.tar.gz` format.

> [!NOTE]
> Since different DB Search and de novo search engines provide different scoring metrics, combining formats within a single Metapepview table may result in skewed or biased results. For this reason, MetaPepView will only add samples to the table if all formats are consistent across samples. However, it is allowed to omit certain import fields for some samples, as long as the formats of present datasets are consistent.
#### DB Search
The dashboard expects a file of peptide matches where every row corresponds to a spectrum match. Multiple db search files can be imported into the field. If *merge db search files* is checked, all files will be concatenated and processed as one sample. The sample name is then specified from the *Name* field. If no db search merge is done, each db search file will be processed as a separate sample. This can be useful if microbial samples will be compared that share the same metagenome (same taxonomy and functional map)

- Peaks 11:  db search spectral match file exported from Peaks Studio 11. The file is usually named `db.psms.csv`.
- Peaks 10:  db search spectral match file exported from Peaks Studio 10. The file is usually named `DB search psm.csv`. Note that Peaks 10 and 11 are incompatible with each other.
- MaxQuant: evidence text file exported from MaxQuant. The file is usually named `evidence.txt`
- Sage (SearchGUI): in SearchGUI, db search analysis can be performed with [Sage](https://github.com/lazear/sage). This results in a file `xxx.sage.tsv`, where `xxx` corresponds to the experiment name. Note that the file is compressed with gzip. Decompress this file first since gzip only compression is not supported by the dashboard (see above note).
#### De Novo
The dashboard expects a file of peptide matches where every row corresponds to a spectrum match. Multiple de novo files can be supplied into the field.

> [!NOTE]
> if *merge db search files* is checked, all de novo files will be processed into a single sample. If it is unchecked (each db search file is it's own sample, see 'DB Search' section), de novo files will be matched to corresponding db search files based on their common MS run. Unmatched de novo files will be ignored.

- Peaks 11:  de novo spectral identification file exported from Peaks Studio 11. The file is usually named `xxx.denovo.csv`, where `xxx` corresponds to the experiment name. 
- Peaks 10:  de novo spectral identification file exported from Peaks Studio 10. The file is usually named `de novo peptides.csv`. Note that Peaks 10 and 11 are incompatible with each other.
- Novor (SearchGUI): de novo identification file generated from Novor in SearchGUI. The file is usually named `xxx.novor.csv`, where `xxx` corresponds to the experiment name. Note that the file is compressed with gzip. Decompress this file first since gzip only compression is not supported by the dashboard (see above note). Also, Novor (SearchGUI) is not compatible with the output generated from [novor.cloud](https://novor.cloud/). This format is currently not supported.
#### Taxonomy Map
The dashboard expects a flat table format with at least a protein name column and a taxonomy id column:
```
Protein_1,tax_a,other_info
Protein_2,tax_b,other_info
Protein_3,tax_b,other_info
Protein_4,tax_a,other_info
```
Protein and taxonomy columns can be defined in the *options* menu. Protein names should be unique (e.g. LCA performed prior) and match to the protein names stored in the db search dataset. One method t

To perform de novo taxonomy annotation of peptide sequences from both db search matched peptides and de novo matched peptides using [Unipept](https://unipept.ugent.be/) API, check the option *global annotation of peptides* in the *options* menu.
#### Functional Map
Protein function datasets that map protein name to KEGG Orthology group id's can be provided here. The dashboard supports two standard formats:

- EggNOG: An annotation file generated from [EggNOG-mapper](http://eggnog-mapper.embl.de/). The file is usually named `xxx.emapper.annotations`.
- gKOALA: tsv file generated from [GhostKOALA](https://www.kegg.jp/ghostkoala/) that maps protein id to KEGG KO. The file is usually named `user_ko.txt`
#### Other options:

- *Export csv*: Export Metapepview table of all samples in csv file format for raw data analysis or manual processing. This dataset cannot be used to import data back into the dashboard. For this, use *Export json*
- *Export json*: Export Metapepview table and metadata as compressed json file format.
- *Import json*: Import Metapepview json back into the dashboard.
- *Experiment name*: Name assigned to export/import files.
#### Import scenario's
Imports are flexibly managed in the dashboard. This means that any import field is only optional. In general, there are two scenario's for import:
##### Import db search data
A sample can be imported from db search data when at least one of the following is included in the import:
- Taxonomy map file
- Functional map file
- global annotation of peptides (taxonomy map options)

de novo data is optional.
##### Import only de novo search data
If no db search data is imported, but only de novo data. A sample can still be created if *global annotation of peptides* is checked. Then, taxonomic classification of peptides with Unipept is performed. This is especially useful for very quick taxonomic composition processing when no protein fasta for db search matching is on hand.
### *Quality Control*: Investigate and benchmark MS spectral data quality

This module contains figures to investigate various spectral parameters from supplied MS experiment data. On the top, A single MS experiment can be imported, with a spectral file in mzML format (Recommended to compress into `.zip` or `.tar.gz` prior to import), as well as a corresponding db search and de novo file from a format of choice.
#### Reference Benchmark
This tab contains figures to benchmark the imported sample against a dataset of experiments. The reference dataset needs to be built first (Not yet published in beta) and can then be imported as json into the dashboard.

### *Taxonomies*: Perform taxonomy composition analysis across samples
This module specifically compares taxonomic composition across samples. The following features are present in this module:
- Taxonomy can be visualized as a stacked barplot or as a heatmap.
- Different taxonomy ranks can be visualized.
- *Normalize abundances*: Represent taxonomy composition of each sample as the fraction of total microbiome composition.
- *Allow global annotation*: for samples that lack protein-taxonomy map import, Unipept based de novo taxonomy classification will be used. Here, both db search data and de novo data is included.
- *Include unannotated*: Add peptides that have no classification at the specified taxonomy rank into the composition into a separate group "Undefined".
- *Filter by Clade*: Show only taxonomy groups that belong to a user specified taxonomy clade.
- *Export taxonomy*: Export taxonomy composition of samples into an easy to process csv. With sample name, taxonomy name, taxonomy rank, peak intensity (if present in data), peptide match count columns. 

### *Taxonomies DB search vs de novo*: Compare taxonomic composition within sample between db search classification and de novo classification
This module tests the coverage and quality of the protein database used for db search matching by comparing the resulting taxonomy composition against a de novo (Unipept) classified taxonomy composition.

- Taxonomy can be visualized as a stacked barplot or as a heatmap.
- Different taxonomy ranks can be visualized.
- *De novo only global annotation*: in de novo (global) taxonomy composition group, omit peptide sequences identified in db search data from the composition processing, considering only de novo identified peptides.
- *Normalize abundances*: Represent taxonomy composition of each sample as the fraction of total microbiome composition.
- *Include unannotated*: Add peptides that have no classification at the specified taxonomy rank into the composition into a separate group "Undefined".
- *Filter by Clade*: Show only taxonomy groups that belong to a user specified taxonomy clade.
- *Export taxonomy*: Export taxonomy composition of samples into an easy to process csv. With sample name, taxonomy name, taxonomy rank, peak intensity (if present in data), peptide match count columns. 
### *Functions*: Visualize protein expression across samples.
This module visualizes abundance of protein functions across samples. Protein functions can represented in various KEGG formats and filtered by [pathway](https://www.genome.jp/kegg/pathway.html) or [BRITE](https://www.genome.jp/kegg/brite.html) groups of interest. In addition, expression data can be exported into a kegg pathway visualization (for pathway groups that support it).

- *KEGG ID display format*: group and quantify proteins based on different id formats (KO, EC, Module, Protein/Gene name).
- *Select groups*: select to filter protein id's to visualize by pathway (and module) group, BRITE group, or by manually selecting KO id's.
- *Normalize sample*: rescale protein function abundance across samples by setting the abundance of a specified sample to 1. All other samples will be displayed as a fold change quantity over the selected normalization sample.
- *Combine multiple annotations*: some peptides may classify to multiple KEGG id's. If checked, the combination of id's is considered a separate group in the figure. If unchecked, multiple annotations will be split and counted towards each group (no peptide will be counted towards any single group twice).
- *Include taxonomies* (Experimental): Show for each sample what taxa contribute towards each protein function expression. (Unstable at the moment).
- *Fractional abundances*: Quantify protein function expression abundance for each sample as fraction compared to total protein quantity.
- *Filter by Clade*: Quantify only protein function expression from peptides that belong to a user specified taxonomy clade.
- *Visualize pathway map*: Given a selected KEGG pathway, select samples in the dropdown menu and by clicking *Show pathway map*, a KEGG page will be loaded showing the pathway with protein functions colored according to the corresponding sample. A color table matching sample to color will shown in the dashboard.
# Referencing

Preprint for the dashboard is currently in preparation.

When using the database for publication, in addition to referencing the tools used as input for the dashboard, the following methods and tools are used during processing

**For de novo (or global) taxonomic annotation of peptide sequences**
- Robbert Gurdeep Singh, Alessandro Tanca, Antonio Palomba, Felix Van der Jeugt, Pieter Verschaffelt, Sergio Uzzau, Lennart Martens, Peter Dawyndt, and Bart Mesuere, "Unipept 4.0: functional analysis of metaproteome data", Journal of Proteome Research, 2019, https://doi.org/10.1021/acs.jproteome.8b00716 
- Kleikamp, H. B., Pronk, M., Tugui, C., da Silva, L. G., Abbas, B., Lin, Y. M., van Loosdrecht, M. C. M. & Pabst, M. (2021). Database-independent de novo metaproteomics of complex microbial communities. Cell Systems, 12(5), 375-383, https://doi.org/10.1016/j.cels.2021.04.003
- Hugo B C Kleikamp, Ramon van der Zwaan, Ramon van Valderen, Jitske M van Ede, Mario Pronk, Pim Schaasberg, Maximilienne T Allaart, Mark C M van Loosdrecht, Martin Pabst, NovoLign: metaproteomics by sequence alignment, ISME Communications, Volume 4, Issue 1, January 2024, ycae121, https://doi.org/10.1093/ismeco/ycae121
