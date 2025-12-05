## Technical overview annotation pipeline

This section provides a high level overview of the processing steps performed during import and annotation of metaproteomics data. Here, details regarding aggregation and quantification of metaproteomics data is outlined. Several decisions were made in missing data handling and format restrictions to manage potential identification and quantification biases from the datasets. The sections outlined below loosely correspond to functions executed in the dashboard backend. Specific processing steps executed during import depends on user specified annotation options and sources that are included into the dashboard (Merge DB search files into one sample, perform Unipept annotation of peptide sequences, etc.).

### Annotation start: Load dataset and set annotation options
1. Check presence of datasets, are datasets present to do annotation? (see scenarios)
2. Load current sample table, Check and filter DB search data in input dataset and current sample table, also, if sample names have duplicates, add suffix to the new import to separate them.
3. Perform [[annotation-pipeline#Annotate peptides|peptide Annotation]] to process data and expand the sample table.
4. Concatenate `MetaPepTable` object for new sample with the current sample table.
### Annotate peptides
1. Import taxonomy db (if GTDB chosen, do not perform Unipept search)
2. Import taxonomy map file (if uploaded) (perform LCA if protein ID present multiple times)
3. Import function map file (if uploaded)
4. Load cRAP dataset if specified by user
5. Import *de novo* data file(s) into MetaPepDeNovo format
	1. Create a dictionary that maps: `{raw_spectral_name: MetaPepDeNovo_obj}`
	2. When importing and processing*de novo* filter out crap peptides if specified
	3. Add *de novo* data only if its `raw_spectral_name` is not present in the dict. Thus, ensure that every spectrum file has only one *de novo* file within a single sample. (between samples, duplicate files are allowed)
6. Import database search file(s) **if present**
	1. If **Merge DB search** is True (all DB search files are one sample):
		1. Load all DB search files to `MetaPepDbSearch` ([[annotation-pipeline#Convert source specific metaproteomics format into meta-PepView format|Load metapep DB search]]) format and store in list, filter out crap peptides if specified for all files
		2.  Concatenate `MetaPepDbSearch` tables (one per DB search file)
		3. Check during concatenation that source files (raw spectrum file) are unique across objects, otherwise, peptides may be counted double
		4. [[annotation-pipeline#Process data into new sample|Create new sample dataset]] for concatenated `MetaPepDbSearch` file.
	2. If **Merge DB search** is False (each DB search file is its own sample):
		1. loop through DB search files:
			1. Load an DB search file and process into `MetaPepDbSearch` (see [[annotation-pipeline#Convert source specific metaproteomics format into meta-PepView format|Load metapep DB search]]) (filter crap if specified)
			2. [[annotation-pipeline#Process data into new sample|Create new sample dataset]] for single `MetaPepDbSearch` file.
			3. Append `MetaPepTable` into list and iterate to next DB search file.
		2. Concatenate all `MetaPepTable` objects into single object..
7. If no DB search files supplied but only *de novo* file: Build `MetaPepTable` (see `build_metapep_table()` level) for sample with only *de novo* data and taxonomy mapping. Append sample to existing `MetaPepTable` and return it.

### Convert source specific metaproteomics format into meta-PepView format
1. Load correct DB search class based on DB search/*de novo* format (e.g. Sage).
2. Read data in class and add sample name to it.
3. Convert format specific DB search/*de novo* object to `MetaPep{DbSearch | DeNovo}` object.
	1. Rename columns to consistent format (save confidence format as variable).
	2. Extract aa sequence from peptides (equate Leucin-Isoleucine, remove PTM).
	3. Filter cRAP peptides out.
	4. Extract all unique source file names to store in a list in the  `MetaPep{DbSearch | DeNovo}` object.
	5. remove file type suffix from source file, format protein ID delimiter.

### Process data into new sample
1. *If DB search data supplied*:
	1. Apply confidence threshold cutoff and aggregate spectrum matches to peptide sequence groups: Create PSM column that counts number of observations of peptide sequence, sum MS1 precursor signal intensities, store maximum spectrum confidence as peptide sequence confidence. From the highest confidence scan, take the spectrum information as peptide information (e.g. retention time, m/z, ppm, scan number, etc.)
	2. `taxonomic_annotation()` (*if protein-taxonomy map present*): supplement peptide grouped data with taxonomy ID and lineage information. If no mapping file present, add empty columns instead. If multiple proteins matched against peptide, store LCA of protein taxa.
	3. `functional_annotation()` (*if protein-function map present*): supplement peptide data with KEGG KO information from function mapper. If multiple proteins mapped against peptide: either, only store information if no conflict in ID between proteins (empty values are ignored), or concatenate IDs into a combined string.
2. *If de novo data supplied*:
	4. `include_de_novo()`: 
		1. *If all DB search files are one single sample*, take all *de novo* file data and concatenate to single `MetaPepDeNovo` (all *de novo* peptides are included in the sample, no matter the spectrum file source)
		2. *If each DB search is separate sample*: fetch *de novo* files that match to source files in DB search data files, concatenate only these files (one DB search file is one sample, this file will only match *de novo* files that come from the same MS runs as that DB search file)
		3. `process_de_novo_data()`: Apply confidence cutoff, peptide length cutoff and group spectra in the concatenated *de novo* object by peptide sequence and sample name. (see `process_db_search_data()` for aggregation rules)
		4. `merge_de_novo_db_search()`: add *de novo* fields to DB search peptide dataset, match by peptide sequence.
		5. Store *de novo* metadata (confidence format, import status, *de novo* file format)
3. *If Unipept taxonomy selected*
	1. `global_taxonomic_annotation()`:
4. Set metadata fields: formats, what data imported, etc.
5. Combine peptide dataset with metadata into `MetaPepTable` object
