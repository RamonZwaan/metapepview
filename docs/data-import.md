# Sample data import

A *project* consists of a dataset of samples, where each sample contains metaproteomics data (peptide sequences from DB search and *de novo* identification), as well as taxonomy and functional information for each peptide. Generally, a "sample" relates to the input sample (metaproteome) used for MS analysis and it may contain peptide data from multiple MS runs (in case of replicates or fractions).

To extend an existing project or start a new one, new samples can be imported into the MetaPepView dashboard. This is done in the [[getting-started#Getting started|Create project]] page. Here, new samples are added from the importer module:

![[Importer_module.PNG]]

For each data component, MetaPepView supports input formats from various sources. The [[prepare-input-data|Prepare input data]] section contains a detailed overview of the input formats that is expected by MetaPepView.

## :warning:{ .yellow } IMPORTANT: Consider the input format constraints

When combining multiple samples into a single sample table, MetaPepView is strict in managing import formats that is allowed. When data is imported into an existing project, MetaPepView performs the following checks if a sample is compatible for import:

- All samples should be compatible in taxonomy and function characterization. Thus, a project can only contain samples with either NCBI or GTDB format taxonomies, not both.
- For each data component, all samples in a project must come from the same data format (e.g. Sage, Novor, GhostKOALA). The output of a sample may differ significantly based on the used input source, to ensure consistency across samples, MetaPepView allows no mixing of input sources. The input formats that the project follows and allows for new data can be seen on the top bar of the *Create project* page.

## Input data components

A sample consists of the following import components:

- DB Search metaproteomics datasets
- *De novo* metaproteomics datasets
- A taxonomy annotation dataset (and/or taxonomy annotation by sequence matching in Unipept)
- A functional annotation dataset

## Import multiple DB search/*de novo* files at once

For DB search and *de novo* metaproteomics input data, multiple files may be selected for import. However, only one taxonomy annotation and functional annotation file may be imported per sample. When multiple DB search or *de novo* metaproteomics files are selected, MetaPepView will either merge the datasets into a single sample, or generate a separate sample for each DB search file. This option is specified in the *merge DB search files* toggle. Merging multiple DB search and *de novo* files allows for easily combining fractions or replicates from a single experimental sample.

## Allowed combinations of input data components 

MetaPepView allows flexible inclusion or omission of input data, it will process all data supplied by the user and ignore processing tasks if a required data source is absent. However, a new sample should contain at least two import components:

- One metaproteomics dataset (peptides from DB search or *de novo* identification).
- One annotation source (taxonomy annotation or functional annotation).

Some import components may be omitted if no data is present. However, this may limit the visualization tools that are available for the sample. In general, the input data configuration depends on the presence or absence of DB search input data.

### Import sample with DB search input

When one (or multiple) DB search input datasets are provided, a sample may be imported if at least one of the taxonomy annotation or functional annotation source is provided. The taxonomy annotation source may be a local accession to taxonomy dataset, annotation of peptide sequences by Unipept, or both. One or multiple *de novo* datasets may be supplemented to the sample.

### Import sample without DB search input

With only *de novo* metaproteomics input data in a sample, MetaPepView can still perform taxonomy analysis on the sample. To import a sample with only *de novo* data, simply check the checkbox *Annotate peptides to Unipept*. MetaPepView can currently only perform taxonomy annotation of peptides by sequence matching in Unipept. No local taxonomy annotation file or functional annotation dataset can be added.

### Import multiple samples at once

MetaPepView allows import of a collection of metaproteomics experiments at once. This is done by loading multiple DB search datasets into the *DB search* field and unchecking the *merge DB search files* toggle. Note that this requires all DB search files to be compatible with the same taxonomy and functional database. In addition, sample names are automatically named after the DB search file names. Thus, these should be named accordingly. This option may be useful for time series analysis or separate import of replicates.

When *de novo* metaproteomics data is included together with DB search data. MetaPepView will internally match each *de novo* dataset with the corresponding DB search dataset through their common spectral raw file names. This ensures that they belong to the same mass spectrometry analysis. *De novo* files that do not match to any DB search dataset are ignored.

!!! note
    Matching *de novo* metaproteomics data to DB search data can get complicated when single DB search data files or *de novo* data files contain multiple raw source files (multiple MS runs merged into one file). Currently, MetaPepView will merge all *de novo* files that share at least one raw source file with the DB search file. Though, this may result in redundancy where one *de novo* peptide identification may be added to multiple samples.

## Import components

### DB search

Metaproteomics data from DB search matching of mass spectrometry data is imported in the *DB search* field. This field expects a dataset of peptide sequence identified scans (one row corresponds to a single spectrum match). MetaPepView supports several DB search output formats:

- [Peaks Studio 10/11](https://www.bioinfor.com/peaks-11/)
- [MaxQuant](https://www.maxquant.org/)
- [Sage](https://sage-docs.vercel.app/)

More information about the expected input format can be found [[prepare-input-data#Import DB search matching data|here]].

Either a single DB search file may be uploaded, or multiple files may be supplied. If *merge DB search files* is toggled on, all files are merged and processed as a single sample. If it is toggled off, each DB search file will be processed as a separate sample. 

In the *Options* menu, a match confidence cutoff may be specified (format depends on DB search format). In addition, a protein ID pattern may be specified to extract from the protein ID fields and make it compatible with Protein ID report formats from other sources. Finally, a pre-filtering of peptide sequences may be performed based on presence in the [cRAP](https://www.thegpm.org/crap/) database. 

### De novo

Metaproteomics data from *de novo* peptide identification of DDA mass spectometry data is imported in the *De novo* field. Similarly to DB search import, this field expects a dataset of peptide sequence identified scans (one row corresponds to a single spectrum identification). MetaPepView supports several *de novo* output formats:

- [Peaks Studio 10/11](https://www.bioinfor.com/peaks-11/)
- [Novor (SearchGUI)](https://compomics.github.io/projects/searchgui)
- [Casanovo](https://casanovo.readthedocs.io/en/latest/)

More information about the expected input format can be found [[prepare-input-data#Import de novo identification data|here]].

Either a single *de novo* file may be uploaded, or multiple files may be supplied. If *merge DB search files* is toggled on, all *de novo* files will be used to supplement the sample dataset together with all DB search data. They will be merged and processed as a single sample. If *merge DB search files* is toggled off, then each *de novo* file will be matched to a DB search file, based on shared raw file names (originating from the same MS run), and processed separately with the DB search data to form a separate sample for each DB search file. *De novo* files that do not match to any DB search dataset are not processed. 

!!! info
    In most input formats, a single DB search or *de novo* file may contain multiple raw file sources. For example, the database search engine used for spectrum matching may combine multiple sample fractions in a single DB search file. The same is true for *de novo* data. In MetaPepView, if DB search files are not merged (each file is a separate sample), supplementing *de novo* data to each DB search file sample is done by matching raw file sources with loose criteria. This means that when one shared raw file source is found between a DB search file and a *de novo* file, they are matched. As a result, one DB search file may be matched to multiple *de novo* files if it contains multiple raw file sources.

In the *Options* menu, a identification confidence cutoff may be specified (format depends on *de novo* format). Also, a pre-filtering of peptide sequences may be performed based on presence in the [cRAP](https://www.thegpm.org/crap/) database. 

!!!WARNING
    Setting an appropriate confidence cutoff is essential for *de novo* identification data, but not trivial. One challenge is that there is no reliable strategy for false discovery rate estimation. Also, only a small fraction of mass spectrometry scans contain the full information to identify a full peptide sequence. However, most *de novo* identification engines will attempt to resolve a peptide sequence from incomplete spectral profiles, resulting in a large fraction of incorrect sequences.

### Taxonomy annotation

The *Taxonomy annotation* field is used for supplementing taxonomic information to the peptides provided in the metaproteomics import data. The main source for taxonomic annotation is from a user supplied accession-taxonomy mapping file (often referred to as *local annotation*). MetaPepView supports import of [GhostKOALA](https://www.kegg.jp/ghostkoala/) format output data (for expected format, see [[prepare-input-data#Taxonomy annotation from GhostKOALA|here]]), or a user generated accession-taxonomy map file. 

For user generated taxonomy files, they should contain an *accession* column, which may be *protein IDs* linked to peptide sequences from the DB search dataset or *peptide sequences* directly. These accessions are mapped with a taxonomy ID, which should be stored in a separate column. MetaPepView supports annotation of NCBI or GTDB format taxonomies. Since the dataset matches protein IDs to taxonomy, this source will only annotate DB search data. More information about the expected format of the protein-taxonomy map is found [[prepare-input-data#Import taxonomy annotation data|here]].

MetaPepView provides an extensive menu with options how it should parse the accession-taxonomy map dataset:
![[Tax_annot_options.PNG]]

- **Delimiter**: The user may set the expected delimiter for the tabular data format.
- **Accession type**: Specify if accessions are protein IDs, or peptide sequence. (MetaPepView will automatically parse the correct data from the DB search dataasets.)
- **Accession parser**: Specify how the accession values, either as complete string or by extracting a substring. It should only be used for protein ID accessions. This can be used to match the accession (protein ID) format to those reported in DB search data.
- **Pattern match (regex)**: Set a *regex pattern* to match a pattern group of the accession column (Only if **Accession parser** is set to "Custom regex").
- **Accession/Taxonomy column index**: Set column indices (starting at 0) for accession and taxonomy values. 
- **Taxonomy element format**: Specify if the taxonomy column contains *taxonomy IDs* or *taxonomy names*.
- **To NCBI taxonomy ID (genome IDs only)**: (GTDB ONLY) For a dataset of genome IDs present in the GTDB database, convert genomes to a NCBI taxonomy ID, then process data with NCBI taxonomy format.

MetaPepView can collect taxonomy information by querying the peptide sequences directly to [Unipept](https://unipept.ugent.be/) (is enabled by the *Annotate peptides to Unipept* switch at the top of the sample import box). This allows taxonomy annotation of both DB search data and *de novo* data, and does not require any additional dataset to be imported. In Unipept, peptide sequences are matched against the UniprotKB protein database, and taxonomy information is collected from matching proteins. MetaPepView requests the *Last Common Ancestor* for each peptide sequence. 

The taxonomy information fetched from Unipept is stored separately from the *local* annotated taxonomy information. Therefore, both strategies may be performed simultaneously for a single sample. This is useful when the local taxonomy map is validated in the *Evaluate community composition* module. It may also be used to supplement the local taxonomy annotation data in the *Community composition* module.

!!! note
    Unipept only supports taxonomy annotation in NCBI format, it cannot be used as an extra taxonomy annotation strategy for GTDB format taxonomy.

!!! warning
    It is not recommended to convert GTDB taxonomy to NCBI format. The conversion may result in a loss of taxonomy information, may decrease the resolution of the taxonomy annotation, and in the worst case, cause inaccurate taxonomy classifications. Use this only if the sample needs to be in NCBI format and there is no other way to obtain NCBI format taxonomy annotations.

### Functional annotation

The *Functional annotation* field is used for supplementing protein functional information to the peptides provided in the metaproteomics input data (only DB search). MetaPepView supports processing and visualization of functional information in [KEGG orthology](https://www.kegg.jp/) format. MetaPepView supports import of functional annotation data from two sources:

- [EggNOG](http://eggnog5.embl.de/#/app/home)
- [GhostKOALA](https://www.kegg.jp/ghostkoala/)

More information about expected format of functional annotation data can be found [[prepare-input-data#Functional annotation|here]].

Both tools perform functional annotation of proteins based on a protein sequence database provided as fasta. Therefore, functional information is mapped to protein IDs, which have to be matched to the protein IDs linked to the peptide sequences inside the DB search files. At the moment, functional annotation is only possible for DB search matched peptide data.

For both datasets, only the *KEGG Orthology* (KO) identifier is of interest for the MetaPepView dashboard. However, EggNOG provides annotations for proteins from several databases (Pfam, COG, gene ontology, *etc.*). These are preserved when exporting the project table in CSV format and can thus be analyzed manually.

!!! note 
    GhostKOALA and EggNOG both use their own databases to functionally annotate proteins. As a result, proteins might be differently annotated depending on what tool is being used. If a particularly complex sample is being analyzed, it might be worthwhile to experiment with both tools to observe the effect on the output.

Options for configuration of functional annotation data is limited. Similarly to DB search and taxonomy import, a protein ID pattern may be specified to extract from the protein ID fields in the functional annotation dataset to make it compatible with Protein ID report formats from DB search data. 

In addition, there is the option, *Combine multiple annotations*, that may be checked. When this option is checked, it will combine functional annotations for peptides that have multiple protein matches: For each field it will get all variants and concatenate them into a new value. When unchecked, it will discard functional annotations of peptides if a conflict is found.
