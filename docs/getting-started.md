# Getting started
When the MetaPepView server is running, the interface is accessed through the web browser from the URL: `http://localhost:8050`. When starting a new session, the dashboard will display the *Create project* page.

![[Startup_screen.PNG]]
*Dashboard interface at startup*

The dashboard contains several modules for data management and visualization. These are accessed in the sidebar. The dashboard contains the following modules:


- *Create project*: Manage project data by importing samples from metaproteomics experiments.
- **Community analysis**:
    - *Community composition*: Visualize community taxonomy composition across samples from the project table.
    - *Community functions*: Visualize functional expression quantification across samples from the project table.
- **Experiment validation**:
    - *Experimental performance*: Visualization toolbox for evaluation of metaproteomics experiment performance metrics, as well as benchmarking of experiments against an experiment reference dataset.
    - *Evaluate community composition*: Compare community taxonomy composition from one sample between peptide annotation from local protein database matching and peptide sequence matching against Uniprot TrEMBL through the Unipept API.

## Download public databases

For import of metaproteomics data and visualization of taxonomy and function groups, MetaPepView requires access to public taxonomy and function databases. These have to be downloaded from their sources. To infer lineages from taxonomy annotations and parse the taxonomy tree, the NCBI and GTDB taxonomy datasets are required (depending on the selected taxonomy annotation format). For visualization of function groups, A dataset mapping KEGG functional groups need to be built.

The presence status of public databases is shown on the bottom of the sidebar. From there, a menu can be accessed to automatically download the required datasets for processing and visualization in the dashboard.

![[download_public_dbs.PNG]]
*Public database download menu*

Once the public databases are imported, a new project table can be created by [[data-import|Importing new samples]]. these samples can be visualized in the taxonomy and functional visualization modules.

![[Import_data.PNG]]
*Build project table by importing metaproteomics experiment data*
