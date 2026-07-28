# Troubleshooting

Meta-PepView relies on data output from many sources, as well as experimental data from a wide variety of experimental setups. As a result, bugs may occur during operation of the dashboard. Most often, this is due to slightly different format from dataset imports compared to what was expected from the dashboard. Most of the time, meta-PepView will provide an informative message what went wrong when encountering an error. If a bug is encountered, please notify by creating an issue on GitHub. Alternatively, send an email to one of the maintainers mentioned on the project page on GitHub.

When sending a bug report, please include all error messages in the dashboard, as well as on the console, if present. Also, provide information about what data files were imported and processed, what sources they come from, and ideally a dataset snippet showing the format of the dataset.

Before starting with meta-PepView, please read the most common problems you may encounter when working with the dashboard below.

## Common issues

### Trying to import a dataset but no response at all from the dashboard

Most likely, a too large dataset was imported for the web browser to accept. In such case, the browser simply ignores the request and does not respond. New versions of the dashboard employ a chunked upload for datasets expected to exceed upload limits (mzML data, project files, etc.), thus this issue should not be present. However, for some datasets, this limit may persist.

If you encounter this issue, first check that you are on the latest version of the dashboard. If the problem persist, send an issue providing information about the dataset you tried to import and the dashboard behavior.

### Trying to import a dataset, but receiving "invalid format" message

Meta-PepView relies on strict data formatting to ensure that all data is correctly parsed. However, source tools that generate output data may change the formatting of the report dataset. Since meta-PepView supports a wide range of data sources, changes in data reporting that breaks compatibility may occur occasionally. If a data import error is encountered, please create an issue on GitHub, providing information about the data source that causes issues, the version of the tool used to create the dataset as well as relevant configuration settings for that tool, and (a snippet of) the dataset, including headers if present.

### Importing sample into the dashboard with "Add sample", but error encountered or no sample added to the project data

Failing to process a sample after successfull import into the dashboard may be due to several potential issues:

- **No valid peptides**: during sample processing, meta-PepView performs several filtering steps. This may potentially result in zero valid sequences. This is usually caused by import setting configuration. In this case, check import settings for unusual values: most notably, ensure that peptide confidence threshold values in DB search and *de novo* data are sensible, as confidence score ranges vary between data source.
- **Inconsistent source formats across samples**: to ensure that meta-PepView compares samples only when they are consistently processed, it forces all samples to follow the same source formats. For example, you may not add a Sage processed sample to a project with MaxQuant samples. When you add a sample with source conflicts, meta-PepView will not add the sample to the project. Ensure that samples you add are consistent with formats outlined in the project data section.
- **Data formatting problems**: Since meta-PepView validates every file imported into the dashboard, the data should be correctly formatted for processing. Nevertheless, some formatting errors may be missed, resulting in processing failure. In that case, please report an issue on GitHub with error information, and information regarding the datasets used (format, tools used, tool versions).

### Dataset imported, but no taxonomy / function profile

If samples have been successfully imported, they will be shown in the "Project data"
table. This table shows datasets that have been imported for each specific sample. Meta-PepView will automatically visualize the samples in taxonomy / function profiling graphs if the required visualization data is present. Otherwise, no information is plotted. If unexpectedly no taxonomy or function information is shown. There may be several causes:

- **The specific data field is not present for the sample**: Some visualization settings requires presence of *de novo* peptides, either local and/or Unipept taxonomy, etc. In other cases, the imported data format does not support some data fields. Notably, signal intensities are not always reported in DB search PSM or *de novo* reports. In this case, look carefully into the Project data table to see what data is imported. You may export the dataset as csv to look further into the raw data.

- **Failed mapping of taxonomy/function datasets and DB search dataset**: DB search peptides are annotated by mapping protein ID's to the taxonomy and function dataset. For any reason, the Protein ID mappings may not correctly correspond between the datasets, for example due to protein ID processing in one of the tools. To check for this, look into the datasets and compare protein ID's between the tool. In case of differences due to ID formatting, you may correct for this by setting regular expression settings in the options pand for the corresponding dataset import.

- **For function profile, KEGG dataset not imported:** The functional annotation profile graph relies on KEGG mapping datasets to process KEGG KO's. If these are not present (shown at the bottom of the side bar), no functional profile figure may be displayed. Meta-PepView should notify the user that no KEGG dataset is present.

### Performance evaluation lacks some performance parameters

The full performance evaluation toolbox relies on a spectral dataset (mzML), feature dataset (featureXML), DB search peptides, and *de novo* peptides to provide all performance metrics. If certain data is missing from the import, then performance metrics that rely on these are not shown. In rare cases, unexpected data formatting may result in absence of performance metrics where they should be present. When encountering absence of performance metrics, first, check if all datasets are imported into the dashboard in the Project data section. If that is the case and data is still missing. Create an issue, providing information about the missing performance metrics, datasets imported, as well as the versions of tools to create the datasets, and a snippet of the data showing the data formatting.

### No Dashboard error, but console shows error/warning messages

Even though the dashboard may operate without problems, the console may display error messages or warnings. This is often due to data processing libraries like Pandas, encountering potential issues in processing. While these may not impact dashboard functionality, it is still important to manage these problems. If such errors are present, please create an issue providing the error information, as well as what dashboard operations induce the error.