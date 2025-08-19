# Processing and annotation

After import data [[data-import|is provided]] into the dashboard, import of the sample into the project table is started by pressing the *Add sample* button. This may be an resource intensive task depending on the size of the input data and may take a minute to execute. During import, MetaPepView performs the following processing tasks:

- Wrangling the input data into a consistent format that is workable for MetaPepView to visualize.
- Grouping and quantifying peptides from scan counts and signal intensities. Grouping is done per sample.
- Checking compatibility of the new sample import sources with those present in the project table.
- Annotating DB search and de novo metaproteomics data with taxonomy information and functional information from the supplied datasets.
- Update the project table with the new sample data.

A high level overview of the importer and annotation pipeline executed during sample import is provided [[annotation-pipeline|here]]. 