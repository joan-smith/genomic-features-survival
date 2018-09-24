# genomic-features-survival
Scripts supporting identification of genomic features affecting survival time in cancer.

## Getting Started

The easiest way to get started is to use `requirements.txt` to set up a [conda](https://conda.io/docs/) environment with the relevant packages.

```$ conda create --name <env> --file requirements.txt```

The code in this repo was built with python 2.7, pandas 0.18, numpy 1.10, and scipy 1.0.0 (as is captured in `requirements.txt`)

The source data is stored in a public [GCS](https://cloud.google.com/storage/) bucket. Documentation for accessing public GCS data is [here](https://cloud.google.com/storage/docs/access-public-data). The bucket for this project is called `public-smith-sheltzer-cancer-analysis`.

Use run.py to download the relevant data from the public GCS bucket *and* perform univariate analyses for copy number and mutation.
run.py takes two optional arguments -- the folder to store data and analysis, (default `.`) and the number of parallel workers to use for the analysis (default none, and a sequential only analysis.). With `-p 4`, run.py takes ~12 hours on a 2017 MacBook Pro.

```$ python run.py -p 4 -o $ouput_directory```


## Organization

 * cbioportal - scripts used to analyze cbioportal data
 * cnv-and-mutations - scripts for analyzing with cnas and and mutations together
 * common-case-zscores - allows getting zscores for every row in a "common" file. common files have genes in rows, patients in columns.
 * common - the common set of utilities and tools used in analysis.
    - `analysis.py` has the meat of cox analysis
    - `mutation_base.py` has the repeatable processing required to turn raw mutation data into usable dataframes.
 * copy-number-analysis - given a copy number file, data about gene/location, and TCGA clinical data, calculate zscores for copy number genes
    - `process_copy_numbers_to_genes.py` has the repeatable processing to turn copy number raw data into usable dataframes.
 * data-munging - contains scripts for miscellaneous small processing tasks: one-off zscores, density plot generation, etc
 * fdr - scripts for performing false discovery correction
 * geo - scripts for analyzing zscores for GEO files
 * make-pancan - given a filetype/platform, take all the per-cancer-type zscore files and produce a file with genes in rows, and cancer types in columns
 * mutation-analysis - given a tcga clinical file, and mutation data from the same set of
patients, calculate zscores and kaplan meier curves for genes mutated in a sufficient number of patients
 * pan-platform - scripts for creating panplatform TCGA files
