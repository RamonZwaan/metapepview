Meta-PepView provides a few pre-built reference datasets to [[experiment-evaluation#Benchmarking|benchmark]] experiments against. However, it is recommended to create an in-house reference database consisting of comparable metaproteomics experiments i.e., samples obtained from similar types of biomaterial, as well as consistent analytical instruments. Meta-PepView provides a command line utility *mpv-buildref*, that allows construction of custom reference datasets for import into the dashboard. This utility extracts metrics from the spectral and (meta)proteomics datasets to provide a consise dataset for the complete range of experiments.

*mpv-buildref* only requires a single root directory location from the user, as well as a few options to manage formatting. *mpv-buildref* will parse all files inside the provided root directory, collect all relevant datasets, and combine all datasets into their respective experiments.

## Running *mpv-buildref*

The *mpv-buildref* tool is available both in the [[installation#Run meta-PepView inside a Docker container|Docker]] version, as well as when installed with [[installation#Install meta-PepView with `pipx` (or `pip`)|`pip`/`pipx`]]. For ease of use, it is recommended to run it from a `pipx` installation, as it will be available directly to the command line:

```Bash
$ mpv-buildref -d peaks11 -n peaks11 -o output.json *root_dir*
```

When running inside Docker, the root directory containing the experiments needs to be mounted to the container file system. Then, the container mount point can be called as the root container when calling *mpv-buildref*. In addition, the output location should be inside a directory (in the container), that is mounted to a host direcory. This may be a separate mount point, or the output location may be inside the root container provided as input:

*After [[installation#Run in Docker|Creating the image]]*:
```Bash
$ docker run -it -v "/host/path/to/root_dir:/home/root_dir" metapepview mpv-buildref -d peaks11 -n peaks11 -o /home/root_dir/output.json */home/root_dir*
```

To get a description of the tool and its options, run it with the `--help` option:

```Bash
$ mpv-buildref --help

usage: buildref [-h] -d {peaks10,peaks11,maxquant,sage} 
                -n {peaks10,peaks11,novor,casanovo} [-o OUTPUT] 
                [-D DB_SEARCH_THRESHOLDS] [-N DE_NOVO_THRESHOLDS] 
                [-I INTENSITY_PERCENTILES] [-T TRANSMISSION_LOSS_PERCENTILES]
                [-r DB_SEARCH_FILE_PATTERN] [-x DE_NOVO_FILE_PATTERN] directory

Create benchmark dataset of metaproteomics experiments from set of experimental 
data. It parses a supplied root directory for all relevant experimental data 
files (mzML, peptide identification data from db search and de novo). All files 
are matched by their common source file name and performance metrics are 
extracted from the data to construct a performance statistics dataset.

positional arguments:
  directory

options:
  -h, --help            show this help message and exit
  -d, --db-search-format {peaks10,peaks11,maxquant,sage}
  -n, --de-novo-format {peaks10,peaks11,novor,casanovo}
  -o, --output OUTPUT
  -D, --db-search-thresholds DB_SEARCH_THRESHOLDS
                        Threshold values used when splitting db search matches 
                        into confidence brackets.
  -N, --de-novo-thresholds DE_NOVO_THRESHOLDS
                        Threshold values used when splitting de novo 
                        identifications into confidence brackets.
  -I, --intensity-percentiles INTENSITY_PERCENTILES
                        MS signal intensity percentiles to be calculated.
  -T, --transmission-loss-percentiles TRANSMISSION_LOSS_PERCENTILES
                        Ion transmission loss percentiles to be calculated.
  -r, --db-search-file-pattern DB_SEARCH_FILE_PATTERN
                        Custom regex pattern used to fetch db search match files 
                        by file name.
  -x, --de-novo-file-pattern DE_NOVO_FILE_PATTERN
                        Custom regex pattern used to fetch de novo 
                        identification files by file name.
```

For correct parsing of proteomics datasets, it is recommended to take the following into account:

- To build a complete dataset, provide for all experiments a spectral file (`mzML`), featuremap (`featureXML`), as well as DB search and *de novo* results.

- When converting raw file names to *mzML* format. **Do not change the name of the mzML file**, and ensure that the *mzML* file name is identical to the raw file name (minus the file type suffix). For each experiment, all data files are automatically linked to each other based on the name of the source file. Depending on the software used for DB search/*de novo* analysis, this will either be the raw spectral data, or the name of the *mzML*. In this case, meta-PepView will assume that the name of the raw file and *mzML* file is the same.

- Ensure that all DB search and *de novo* datasets can be recognized as such from the file names (for example, adding a fixed prefix / suffix to each file). If all experiment results were exported using default filenames provided by the analysis software, the tool should be able to automatically recognize all files inside the supplied root directory. However, if custom file names were assigned, regular expression patterns may need to be supplied for the tool to recognize the datasets (DB search: `-r`, *de novo*: `-x`). If no patterns are supplied, it will search DB search and *de novo* datasets using the following patterns:
    
    **DB search:**

    - peaks11:   *"db.psms.csv"*
    - peaks10:   *"DB search psm.csv"*
    - maxquant:   *"evidence.txt"*
    - sage':       *"\*.sage.tsv"*

    **de novo**

    - peaks11:   *"\*.denovo.csv"*
    - peaks10:   *"de novo peptides.csv"*
    - novor:      *"\*.novor.csv"*
    - casanovo:   *"\*.mztab"*

- It is recommended to set reasonable confidence thresholds for the used DB search (option: `-D`) and *de novo* (option: `-N`) output format. If no thresholds are given, it will default to the following values: DB search: "30 50 80", *de novo*: "50 80 90". While many proteomics tools score identifications in a 0 - 100 range, some use a different range. For example, Casanovo scores matches between -1 - 1. For such cases, custom threshold values should be provided to partition peptides in confidence groups.
