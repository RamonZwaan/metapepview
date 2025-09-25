# Getting started
When the MetaPepView server is running, the interface is accessed through the web browser from the URL: `http://localhost:8050`. When starting a new session, the dashboard will display the *Create project* page.

![[Startup_screen.PNG]]
*Dashboard interface at startup*

The dashboard contains several modules for data management and visualization. These are accessed in the sidebar. The dashboard contains the following modules:


- *Create project*: Manage project data by importing samples from metaproteomics experiments.
- **Community analysis**:
    - *Community composition*: Visualize community taxonomy composition across samples from the project table.
    - *Community functions*: Visualize functional expression quantification across samples from the project table.
- **Experiment evaluation**:
    - *Experiment performance evaluation*: Visualization toolbox for evaluation of metaproteomics experiment performance metrics, as well as benchmarking of experiments against an experiment reference dataset.
    - *Community composition evaluation*: Compare community taxonomy composition from one sample between peptide annotation from local protein database matching and peptide sequence matching against Uniprot TrEMBL through the Unipept API.

## Download public databases

For import of metaproteomics data and visualization of taxonomy and function groups, MetaPepView requires access to public taxonomy and function databases. These have to be downloaded from their sources. To infer lineages from taxonomy annotations and parse the taxonomy tree, the NCBI or GTDB taxonomy datasets are required (depending on the selected taxonomy annotation format). For visualization of function groups, the KEGG function mapping dataset has to be fetched online.

The presence status of public databases is shown on the bottom of the sidebar. From there, a menu can be accessed to automatically download the required datasets for processing and visualization in the dashboard.

!!! note
    When downloading datasets, MetaPepView will create a directory `.metapepview` inside the user home directory. The directory path that MetaPepView will look for and write datasets can be changed on the top of the download menu. Note that MetaPepView will only create directories for each dataset when the user checks the *Create parent dir* boxes. Make sure to check these boxes when downloading the datasets for the first time.

![[download_public_dbs.PNG]]
*Public database download menu*

!!! note
    If MetaPepView is running inside a [[installation#Run in Docker.|Docker container]], downloaded databases will only persist during the lifetime of the container. To ensure that databases will be preserved across container instances, the storage directory can be mounted to a directory on the host computer or a [volume](https://docs.docker.com/engine/storage/volumes/) using the `-v` flag when creating a new container:
    ```Bash
    docker run -p 8050:8050 -v "ref_db_volume:/home/.metapepview" metapepview
    ```

Once the public databases are imported, a new project table can be created by [[data-import|Importing new samples]]. these samples can be visualized in the taxonomy and functional visualization modules.

![[Import_data.PNG]]
*Build project table by importing metaproteomics experiment data*
