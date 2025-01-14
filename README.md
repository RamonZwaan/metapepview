
# Description
MetaPepView (Beta) is a metaproteomics dashboard designed for visualizing the taxonomic composition and expressed enzymatic or metabolic functions from microbial community data, as well as evaluating the performance of metaproteomics experiments. This beta version features a user-friendly, web-based interface that can be run locally on desktop PCs. It enables users to compare results from database searches and de novo annotations while evaluating the quality and coverage of reference databases. Feedback and contributions are welcome as we continue to enhance MetaPepView. Contact: r.vanderzwaan@tudelft.nl

Main features of MetaPepView include:

- Visualize taxonomic compositions and functional profiles of metaproteomics experiments as obtained from various database search or de novo sequencing engines.
- Perform de novo taxonomy classification using the Unipept API and compare the obtained de novo taxonomic composition with db search results.
- Analyze mass spectrometric performance from mzML raw files and evalute db search and de novo matches from spectral data.
- Benchmark experimental performance with proteomcis reference datasets.
- NOTE: the beta version of MetaPepView has been tested for mass spectrometric raw data obtained from QE Orbitrap mass spectrometers. Data conversion to mzML was performed using msconvert https://proteowizard.sourceforge.io/tools/msconvert.html.
  
# Installation / Set up
## PC requirements
- The dashboard is cross-platform and has been tested in Windows 10/11 and Linux.
- A minimum RAM of 16 GB is required for full operation of the dashboard, 32 GB is recommended.
- The dashboard can be accessed through both Firefox (Gecko browser engine) and any Chromium based browser (Chrome, Edge, Brave, etc.). **Note: the Chomium browser engine supports file import up to ~300 MB. Any larger file (common in mzML) will not be imported. Firefox does not have such a limit and multi GB files have been successfully imported.**
- Internet connection is required to download taxonomy and function datasets, and for connection to the Unipept API for de novo taxonomy classification.
  
## Setting up the server
MetaPepView runs on top of the dash framework in python. To run the dashboard, the server needs to be started locally on the PC. The simples approach is to start a [Docker](https://www.docker.com/) container from the DockerFile provided in the Github repository. For this, docker (or another OCI compliant manager) has to be installed on the system.

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

Alternatively the repository provides a requirements file to install all required python packages with [pip](https://packaging.python.org/en/latest/tutorials/installing-packages/#requirements-files). To execute the tool directly on the host pc, a recent version of python needs to be installed (at least 3.11). Then, The dashboard can be run directly on the host pc by executing:
```Shell
$ pip install -r ./requirements.txt
$ python ./index.py
```

# Dashboard operation

## Getting started

> [!WARNING]
> The Dash framework is stateless. This means that when you refresh the web page, all data stored in memory will be wiped and a new session will be started. Therefore, make sure to export (save) experiments after raw data have been imported to the dashboard, with *Export json*. By that, saved data can also be reimported quickly in a new session using *Import json*.

When the MetaPepView server is running and you open it in the browser at `http://localhost:8050`, you first will see the *Annotation* page. Here, you can import the samples into MetaPepView. For every sample you can load db search results, de novo peptide sequencing data, and taxonomic + functionional classification data. Multiple samples can be loaded and compared.
### Import online reference datasets.
Mapping taxonomy IDs to phylogenetic lineages and functions to pathways and terms, reference databases need to be loaded to MetaPepView. This can be done **on the left sidebar at the bottom** at *Database loaded*. By clicking *Fetch source*, you can download the latest version of NCBI taxonomy, GTDB taxonomy and KEGG data from their original sources. When using the dashboard for the first time, select the datasets for download and **create parent directories**. Note, thise databases require approx. **1-5 GB** storage memory.

### *Annotation*: Import new sample
The MetaPepView dashboard is very flexible and allows you to import database search, de novo and taxonomic and functional classification data from various sources:
> [!NOTE]
> For efficient RAM handling, the dashboard also supports reading of compressed datasets. For speed improvement and lower RAM usage, compress large files (>100 MB) in `.zip` or `.tar.gz` format.

> [!NOTE]
> Different DB search and de novo search engines employ different scoring metrices. For this reason, MetaPepView will only allow to import samples that use a consistent type of database search, de novo sequencing and taxonomic and functional classification formats. However, it is allowed to omit certain import fields, as long as this is consistent across the compelte dataset.
#### DB Search
Metapepview does not perfrom database searching itself. However, MetaPepView is very flexible, and it only expects a file where every row corresponds to a spectrum match. Results from different database engines can be imported. If *merge db search files* is checked, all files will be concatenated and processed as one sample. The sample name is then specified from the *Name* field. If no db search merge is done, each db search file will be processed as a separate sample. This can be useful if samples are compared that share the same metagenome (same taxonomy and functional map)

- Peaks 11:  db search spectral match file exported from Peaks Studio 11. The file is usually named `db.psms.csv`.
- Peaks 10:  db search spectral match file exported from Peaks Studio 10. The file is usually named `DB search psm.csv`. Note that Peaks 10 and 11 are incompatible with each other.
- MaxQuant: evidence text file exported from MaxQuant. The file is usually named `evidence.txt`
- Sage (SearchGUI): in SearchGUI, db search analysis can be performed with [Sage](https://github.com/lazear/sage). This results in a file `xxx.sage.tsv`, where `xxx` corresponds to the experiment name. Note that the file is compressed with gzip. Decompress this file first since gzip only compression is not supported by the dashboard (see above note).
#### De Novo
Metapepview does also not perform de novo sequeicnt. However, de novo sequencing outputs from different tools can be imported. Metapepview expects a file containing peptide sequences where every row corresponds to a new sequence. Multiple de novo files can be supplied into the field.

> [!NOTE]
> if *merge db search files* is checked, all de novo files will be combined to a single sample. If it is unchecked (each db search file is it's own sample, see 'DB Search' section), de novo files will be matched to corresponding db search files based on the MS raw data file they originate from. Unpaired de novo files will be ignored.

- Peaks 11:  de novo spectral identification file exported from Peaks Studio 11. The file is usually named `xxx.denovo.csv`, where `xxx` corresponds to the experiment name. 
- Peaks 10:  de novo spectral identification file exported from Peaks Studio 10. The file is usually named `de novo peptides.csv`. Note that Peaks 10 and 11 are incompatible with each other.
- Novor (SearchGUI): de novo identifications generated from Novor in SearchGUI. The file is usually named `xxx.novor.csv`, where `xxx` corresponds to the experiment name. Note that the file is compressed with gzip. Decompress this file first since gzip only compression is not supported by the dashboard (see above note). Also, Novor (SearchGUI) is not compatible with the output generated from [novor.cloud](https://novor.cloud/). This format is currently not supported.
  
#### Taxonomy Map
The dashboard expects a flat table format with at least a protein name column and a taxonomy id column:
```
Protein_1,tax_a,other_info
Protein_2,tax_b,other_info
Protein_3,tax_b,other_info
Protein_4,tax_a,other_info
```
Protein and taxonomy columns can be defined in the *options* menu. Protein names should be unique (e.g. LCA performed prior) and match to the protein names stored in the db search dataset. One method t

To perform de novo taxonomy annotation of peptide sequences obtained from db searching and de novo sequencing using the [Unipept](https://unipept.ugent.be/) API, check the option *global annotation of peptides* in the *options* menu.

#### Functional Map
Functional annotation data that link proteins to KEGG Orthologiess can be provided here. The dashboard supports two input formats:

- EggNOG: annotation file generated from [EggNOG-mapper](http://eggnog-mapper.embl.de/). The file is usually named `xxx.emapper.annotations`.
- gKOALA: tsv file generated from [GhostKOALA](https://www.kegg.jp/ghostkoala/) that maps protein id to KEGG KO. The file is usually named `user_ko.txt`
#### Other options:

- *Export csv*: Export Metapepview table of all samples in csv file format for raw data analysis or manual processing. This dataset cannot be used to import data back into the dashboard. For this, use *Export json*
- *Export json*: Export Metapepview table and metadata as compressed json file format.
- *Import json*: Import Metapepview json back into the dashboard.
- *Experiment name*: Name assigned to export/import files.
  
#### Import scenario's
Imports are flexibly managed, meaning that import fields are optional. In general, there are two scenario's for imports:

##### Import db search data
A sample can be imported from db search data when at least one of the following is included in the import:
- Taxonomy map file
- Functional map file
- global annotation of peptides (taxonomy map options)
de novo data are optional.

##### Import only de novo search data
Only de novo sequencing data are provided. To create a sample check *global annotation of peptides*. Then, taxonomic classification of peptides with Unipept is performed. This is useful for performing a rapid taxonomic profiling of the sample.

### *Quality Control*: Investigate and benchmark experimental performance
This module provides figures and various spectral parameters from mass spectrometric mzML raw data and search data. On the top, a single experiment can be loaded, with a spectral file in mzML format (Recommended to compress into `.zip` or `.tar.gz` prior to import), as well as a corresponding db search and de novo sequencing file.

> [!NOTE]
> Quality control is currently designed for DDA (Data-dependent Acquisition) experiments. The quality control module has been tested on mzML files derived from QE Orbitrap raw data. mzML conversion performed with [MSConvert](https://proteowizard.sourceforge.io/index.html).

#### Reference Benchmark
This module benchmarks the imported sample against a set of reference data with varying quality levels (from high to low). However, the reference dataset must be created (not yet included in the beta version) and imported as a JSON file before it becomes available in MetaPepView.

### *Taxonomies*: analyse taxonomic composition across (multiple) samples
This module allows you to analyze and compare the taxonomic composition of one or multiple samples. It offers the following features:
- Taxonomic composition shown as a stacked barplot or as a heatmap.
- Different taxonomic ranks can be visualized and compared.
- *Normalize abundances*: show composition of each sample as fraction of total (100%).
- *Allow global annotation*: where no taxonomic annotations are imported, the de novo taxonomy will be created using Unipept. Here, both db search data and de novo data are considered.
- *Include unannotated*: add peptides that have no classification at the specified taxonomic rank into the composition into a separate group "Undefined".
- *Filter by Clade*: Show only taxonomies that belong to certain taxonomic group, as specified by the user.
- *Export taxonomy*: Export taxonomic composition meta data to csv file. This usually inclues, the sample name, taxon names, taxonomic ranks, peak intensities (if available), peptide match counts. 

### *Taxonomies DB search vs de novo*: Compare taxonomic composition as obtained by db searching and de novo annotation
This module compares the taxonomic composition derived from database searching with that obtained through de novo annotation. It enables to evaluate quality and coverage of the database used for database searching. This module also flags organisms that are absent from the database.

- Taxonomy can be visualized as a stacked barplot or as a heatmap.
- Different taxonomic ranks can be visualized.
- *De novo only global annotation*: in de novo (global) taxonomy composition group, omit peptide sequences identified in db search data from the composition processing, considering only de novo identified peptides.
- *Normalize abundances*: Represent taxonomy composition of each sample as the fraction of total microbial abundance.
- *Include unannotated*: Add peptides that have no classification at the specified taxonomy rank into the composition into a separate group "Undefined".
- *Filter by Clade*: Show only taxonomies that belong to a taxonomic group, as specified by the user.
- *Export taxonomy*: Export taxonomy composition of samples into an easy to process csv. With sample name, taxonomy name, taxonomy rank, peak intensity (if present in data), peptide match count columns.
   
### *Functions*: Visualize protein expression across samples.
This module visualizes enzymatic and metabolic functions. Protein functions can represented in various KEGG formats and filtered by [pathway](https://www.genome.jp/kegg/pathway.html) or [BRITE](https://www.genome.jp/kegg/brite.html) groups of interest. In addition, expression data can be exported into a kegg pathway visualization (for pathway groups that support it).

- *KEGG ID display format*: group and quantify proteins based on different id formats (KO, EC, Module, Protein/Gene name).
- *Select groups*: select to filter protein id's to visualize by pathway (and module) group, BRITE group, or by manually selecting KO id's.
- *Normalize sample*: rescale protein function abundance across samples by setting the abundance of a specified sample to 1. All other samples will be displayed as a fold change quantity over the selected normalization sample.
- *Combine multiple annotations*: some peptides may classify to multiple KEGG id's. If checked, the combination of id's is considered a separate group in the figure. If unchecked, multiple annotations will be split and counted towards each group (no peptide will be counted towards any single group twice).
- *Include taxonomies* (Experimental): Show for each sample what taxa contribute towards each protein function expression. (Unstable at the moment).
- *Fractional abundances*: Quantify protein function expression abundance for each sample as fraction compared to total protein quantity.
- *Filter by Clade*: Quantify only protein functions from peptides that belong to a user specified taxonomic group.
- *Visualize pathway map*: Given a selected KEGG pathway, select samples in the dropdown menu and by clicking *Show pathway map*, a KEGG page will be loaded showing the pathway with protein functions colored according to the corresponding sample. A color table matching sample to color will shown in the dashboard.

# Referencing
In articles, you can cite as follows: Ramon van der Zwaan, Mark van Loosdrecht and Martin Pabst (2025, preprint in preparation), TU Delft, The Netherlands. Contact: r.vanderzwaan@tudelft.nl

In addition the following articles and tools are used during processing:
**For de novo (or global) taxonomic annotation of peptide sequences**
- Gurdeep Singh R, Tanca A, Palomba A, Van der Jeugt F, Verschaffelt P, Uzzau S, Martens L, Dawyndt P, Mesuere B. Unipept 4.0: functional analysis of metaproteome data. Journal of proteome research. 2018 Nov 22;18(2):606-15. https://doi.org/10.1021/acs.jproteome.8b00716 
- Kleikamp HB, Pronk M, Tugui C, da Silva LG, Abbas B, Lin YM, van Loosdrecht MC, Pabst M. Database-independent de novo metaproteomics of complex microbial communities. Cell Systems. 2021 May 19;12(5):375-83. https://doi.org/10.1016/j.cels.2021.04.003
- Kleikamp HB, Van der Zwaan R, Van Valderen R, Van Ede JM, Pronk M, Schaasberg P, Allaart MT, Van Loosdrecht MC, Pabst M. NovoLign: metaproteomics by sequence alignment. ISME communications. 2024 Jan;4(1):ycae121. https://doi.org/10.1093/ismeco/ycae121
