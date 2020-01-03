# Sex Prediction from RNASeq

The 01-clean_split_data.R, 02-train_elasticnet.R, 03-evaluate_model.R, 04-present_results.Rmd pipeline trains and evaluates an elasticnet logistic regression model to predict sex from RNASeq data.  The training features are gene expression transcripts, and the training labels are reported_gender values for each sample.

The pipeline is a response to issue [#84](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/84). Based on the accuracy achieved here, this classifier can be helpful to predict values for datasets without annotated sex information.

See the comparison results plots [here](https://cansavvy.github.io/openpbta-notebook-concept/snv-callers/compare_snv_callers_plots.nb.html).

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [How to run the caller consensus analysis](#how-to-run-the-caller-consensus-analysis)
- [Summary of Methods](#summary-of-methods)
  - [Variant Allele Fraction Calculation](#variant-allele-fraction-calculation)
  - [Mutation Comparisons](#mutation-comparisons)
  - [Tumor Mutation Burden Calculation](#tumor-mutation-burden-calculation)
    - [All mutations TMB](#all-mutations-tmb)
    - [Coding only TMB](#coding-only-tmb)
- [General usage of scripts](#general-usage-of-scripts)
  - [01-setup_db.py](#01-setup_dbpy)
  - [02-merge_callers.R](#02-merge_callersr)
  - [03-calculate_tmb.R](#03-calculate_tmbr)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## How to run the caller consensus analysis

To run the evaluations and comparisons of all the SNV callers, call the bash script:

```
bash run_caller_analysis.sh
```
This bash script will return:

- Comparison plots in a notebook: [`compare_snv_callers_plots.nb.html`](https://cansavvy.github.io/openpbta-notebook-concept/snv-callers/compare_snv_callers_plots.nb.html).
- A zip file containing:
  - `pbta-snv-consensus-mutation.maf.tsv` - is  [MAF-like file](#consensus-mutation-call) that contains the snvs that were called by all three of these callers for a given sample are saved to this file.
  These files combine the [MAF file data](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) from 3 different SNV callers: [Mutect2](https://software.broadinstitute.org/cancer/cga/mutect), [Strelka2](https://github.com/Illumina/strelka), and [Lancet](https://github.com/nygenome/lancet).
  See the methods on the callers' settings [here](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#somatic-single-nucleotide-variant-calling) and see [the methods of this caller analysis and comparison below](#summary-of-methods).  
  - `pbta-snv-consensus-mutation-tmb-coding.tsv` - Tumor Mutation burden calculations using *coding only* mutations use the consensus of Lancet, Mutect2, and Strelka2.
  - `pbta-snv-consensus-mutation-tmb-all.tsv` - Tumor Mutation burden calculations using *all* mutations use the consensus of Mutect2, and Strelka2. (Lancet was excluded because it has a [coding region bias in the way it was run](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#snv-and-indel-calling)).

## Summary of Methods

### Variant Allele Fraction Calculation

Calculate variant allele fraction (VAF) for each variant.
This is done in `01-setup_db.py`.

```
vaf = (t_alt_count) / (t_ref_count + t_alt_count)
```
This is following the [code used in
`maftools`](https://github.com/PoisonAlien/maftools/blob/1d0270e35c2e0f49309eba08b62343ac0db10560/R/plot_vaf.R#L39).

### Mutation Comparisons

The default consensus mutations called are those that are shared among all of Strelka2, Mutect2, and Lancet.
Mutations were considered to be the same if they were identical in the following field: `Chromosome`, `Start_Position`, `Reference_Allele`,  `Allele`, and `Tumor_Sample_Barcode`.
As Strelka2 does not call multinucleotide variants (MNV), but instead calls each component SNV as a separate mutation, MNV calls from Mutect2 and Lancet were separated into consecutive SNVs before comparison with Strelka2.

### Tumor Mutation Burden Calculation

For each experimental strategy and TMB calculation, the intersection of the genomic regions effectively being surveyed are used.
These genomic regions are used for first filtering mutations to these regions and then for using the size in bp of the genomic regions surveyed as the TMB denominator. 

#### All mutations TMB

For all mutation TMBs, Lancet is not used because of the [coding bias in the way it was run.](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#snv-and-indel-calling)
For WGS samples, the size of the genome covered by the intersection of Strelka2 and Mutect2's surveyed areas is used for the denominator.
```
WGS_all_mutations_TMB = (total # snvs called by Strelka2 and Mutect) / intersection_strelka_mutect_genome_size
```
For WXS samples, the size of the genome the WXS bed region file is used for the denominator.
```
WXS_all_mutations_TMB = (total # snvs called by Strelka2 and Mutect2 ) / wxs_genome_size
```
#### Coding only TMB

Coding only TMB uses all three callers: Strelka2, Mutect2, and Lancet and the intersection demoninators are calculated by using coding sequence ranges in the gtf from Gencode 27.
This file is included in the data download.
SNVs outside of these coding sequences are filtered out before being summed and used for TMB calculations like such:

```
WGS_coding_only_TMB = (total # coding sequence snvs called by all three of Strelka, Lancet, and Mutect2 ) / intersection_strelka_lancet_mutect_CDS_genome_size
```
Because the same WXS BED file applies to all callers, that file is intersected with the coding sequences for filtering and for determining the denominator. 
```
WXS_coding_only_TMB = (total # coding sequence snvs called by all three of Strelka, Lancet, and Mutect2 ) /
intersection_wxs_CDS_genome_size
```

## General usage of scripts

**Overall notes about these scripts:**
- The scripts are sequential as noted by their number.
- All file path-related options assume the file path given is relative to `OpenPBTA-analysis/analyses/sex-prediction-from-RNASeq`.
- The scripts create user-specified output directories that do not exist at run time.
- The scripts add files whose names do not match existing files in the output directories.  Output files whose names match existing files overwrite the previous versions.

### run-sex-prediction-from-RNASeq.sh

A bash shell script that runs the entire pipeline.  Global arguments for the entire pipeline are specified in the USER-SPECIFIED ARGUMENTS section of the script


**Argument descriptions**
```
PROCESSED
output directory of script 01, input directory of scripts 02 and 03

MODELS
outtput directory of script 02, input directory of script 03

RESULTS
output directory of script 03

SEED
seed for R's pseudorandom number generator

FILENAME_LEAD
argument for script 01 output file specification and script 02 and 03 input and output file specification

TRANSCRIPT_TAIL_PERCENT_ARRAY
a bash shell script array to control median absolute deviation filtering of training transcripts.  Array values must be between 0 and 1.  
No limit on the size of the array.  For each value in the array, a model is trained and evaluated by filtering out the bottom (1 - TRANSCRIPT_TAIL_PERCENT_ARRAY value) of transcripts from the training set.

TRAIN_PERCENT
a value greater than 0 and less than or equal to 1.  The expression dataset is partitioned TRAIN_PERCENT for training and (1 - TRAIN_PERCENT) for testing.  If TRAIN_PERCENT = 1, no test set is created.
  
TRAIN_TARGET_COLUMN
argument to specify what column in the target data frame will be used as labels during training

targetColumns
a bash shell script array specifying the columns to be used as test labels when evaluating a model's performance.Performance reports are produced for each value in the array.  No limit on the size of the array.

```

### 01-clean_split_data.R

This script cleans the gene expression data -- i.e., it drops anything with invalid labels in either reported_gender or 
germline_sex_estimate, and if the training partition size < 1.0, it splits the data into training and testing sets that 
are saved as separate files to be used downstream.


**Argument descriptions**
```
  --expression ../../data/pbta-gene-expression-kallisto.stranded.rds 
    File path and file name of the OpenPBTA gene expression file used for the current run of the pipeline(RDS)
  --metadata ../../data/pbta-histologies.tsv 
    File path and file name of the OpenPBTA histologies file used for the current run(TSV)
  --output_directory PROCESSED 
    File path where you would like the cleaned training and test sets to be stored
  --train_expression_file_name TRAIN_EXPRESSION_FILE_NAME 
    
  --test_expression_file_name $TEST_EXPRESSION_FILE_NAME \
  --train_targets_file_name $TRAIN_TARGETS_FILE_NAME  \
  --test_targets_file_name $TEST_TARGETS_FILE_NAME \
  --full_targets_file_name $FULL_TARGETS_FILE_NAME \
  --seed $SEED \
  --train_percent $TRAIN_PERCENT

  -d DB_FILE, --db-file DB_FILE
     Path of the database file to use or create. Defaults to `data.sqlite`.
   --strelka-file STRELKA_FILE
     Path of the MAF formatted data file from the strelka2 caller(TSV).
   --mutect-file MUTECT_FILE
     Path of the MAF formatted data file from the mutect2 caller(TSV).
   --lancet-file LANCET_FILE
     Path of the MAF formatted data file from the lancet caller(TSV).
   --vardict-file VARDICT_FILE
     Path of the MAF formatted data file from the vardict caller(TSV).
   --meta-file META_FILE, --hist-file META_FILE
     Path of the metadata/histology data file(TSV).
   --overwrite           Overwrite tables that may already exist.
```

### 02-merge_callers.R

Using the database created by `01-setup_db.py`, merge callers' data files into consensus [MAF-like file](#snv-caller-comparison-analysis).

**Argument descriptions**
```
 --db_file : Path to sqlite database file made from 01-setup_db.py
 --output_file : File path and file name of where you would like the MAF-like
                 output from this script to be stored.
 --vaf_filter: Optional Variant Allele Fraction filter. Specify a number; any
               mutations with a VAF that are NA or below this number will be
               removed from the vaf data.frame before it is saved to a TSV file.
 --overwrite : If TRUE, will overwrite any reports of the same name. Default is
              FALSE
```
### 03-calculate_tmb.R

Using the consensus file created in `02-merge_callers.R`, calculate TMB for all
WGS and WXS samples.
Two TMB files are created, one including *all snv* called by Strelka2 and Mutect2 (Lancet is excluded from this TMB calculation consensus because of a [coding region bias in the way it was ran](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#snv-and-indel-calling)), and a *coding snvs only* TMB calculation.

**Argument descriptions**
```
 --consensus : File path to the MAF-like file.
 --db_file : Path to sqlite database file made from 01-setup_db.py
 --metadata : Relative file path to MAF file to be analyzed. Can be .gz compressed.
              Assumes file path is given from top directory of 'OpenPBTA-analysis'.
 --bed_wgs : File path that specifies the caller-specific BED regions file.
             Assumes from top directory, 'OpenPBTA-analysis'.
 --bed_wxs : File path that specifies the WXS BED regions file. Assumes file path
             is given from top directory of 'OpenPBTA-analysis'
 --overwrite : If specified, will overwrite any files of the same name. Default is FALSE.
```
