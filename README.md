# genomic-features-survival
Scripts supporting identification of genomic features affecting survival time in cancer


 * cbioportal - scripts used to analyze cbioportal data
 * cnv-and-mutations - scripts for analyzing with cnas and and mutations together
 * common-case-zscores - allows getting zscores for every row in a "common" file. common files have genes in rows, patients in columns.
 * common - the common set of utilities and tools used in analysis.
    - analysis.py has the meat of cox analysis
 * copy-number-analysis - given a copy number file, data about gene/location, and TCGA clinical data, calculate zscores for copy number genes
 * data-munging - contains scripts for miscellaneous small processing tasks: 1 off zscores, density plot generation, etc
 * fdr - scripts for performing false discovery correction
 * geo - scripts for analyzing zscores for GEO files
 * make-pancan - given a filetype/platform, take all the per-cancer-type zscore files and produce a file with genes in rows, and cancer types in columns
 * mutation-analysis - given a tcga clinical file, and mutation data from the same set of
patients, calculate zscores and kaplan meier curves for genes mutated in a sufficient number of patients
 * pan-platform - scripts for creating panplatform TCGA files
