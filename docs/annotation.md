# Processing and annotation for the compositions and functions dataset

After [[data-import#Import samples into compositions and functions dataset|import data is provided]] for a new sample in the *compositions and functions dataset*, import of the sample into the dataset is started by pressing the *Add sample* button. This may be an resource intensive task depending on the size of the input data and may take a minute to execute. During import, meta-PepView performs the following processing tasks:

- Wrangling the input data into a consistent format that is workable for meta-PepView to visualize.
- Grouping and quantifying peptides from scan counts and signal intensities. Grouping is done per sample.
- Checking compatibility of the new sample import sources with those present in the project table.
- Annotating DB search and *de novo* metaproteomics data with taxonomy information and functional information from the supplied datasets.
- Update the project table with the new sample data.

A high level overview of the importer and annotation pipeline executed during sample import is provided [[annotation-pipeline|here]]. 