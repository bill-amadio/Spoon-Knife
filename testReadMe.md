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
  - [run-sex-prediction-from-RNASeq.sh](#run-sex-prediction-from-RNASeq)
  - [01-clean_split_data.R](#01-clean-split-data)
  - [02-train_elasticnet.R](#02-train_elasticnet)  
  - [03-evaluate_model.R](#03-evaluate_model)
  - [04-present_results.Rmd](#04-present_results)

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
- The scripts add files whose names do not match existing files, if any, in the output directories.  Output files whose names 
  match existing files overwrite the previous versions.

### run-sex-prediction-from-RNASeq

A bash shell script that runs the entire pipeline.  Global arguments for the pipeline are specified in the USER-SPECIFIED ARGUMENTS section of the script


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
a bash shell script array to control median absolute deviation filtering of training transcripts.  Array 
values must be between 0 and 1.  No limit on the size of the array.  For each value in the array, a model 
is trained and evaluated by filtering out the bottom (1 - TRANSCRIPT_TAIL_PERCENT_ARRAY value) of transcripts 
from the training set.

TRAIN_PERCENT
a value greater than 0 and less than or equal to 1.  The expression dataset is partitioned TRAIN_PERCENT 
for training and (1 - TRAIN_PERCENT) for testing.  If TRAIN_PERCENT = 1, no test set is created.
  
TRAIN_TARGET_COLUMN
argument to specify what column in the target data frame will be used as labels during training

targetColumns
a bash shell script array specifying the columns to be used as ground truth labels when evaluating a model's 
performance.  Performance reports are produced for each value in the array.  No limit on the size of 
the array.

```

### 01-clean-split-data

This script cleans the gene expression data -- i.e., it drops anything with invalid labels in either 
reported_gender or germline_sex_estimate, and if the training partition size < 1.0, it splits the data 
into training and testing sets that are saved as separate files to be used downstream.


**Argument descriptions**
```
  --expression ../../data/pbta-gene-expression-kallisto.stranded.rds 
    File path and file name of the OpenPBTA gene expression file used for the current run of the pipeline(RDS)
  --metadata ../../data/pbta-histologies.tsv 
    File path and file name of the OpenPBTA histologies file used for the current run(TSV)
    
Values for the following arguments are computed in run-sex-prediction-from-RNASeq.sh using values from the 
USER-SPECIFIED-ARGUMENTS.

  --output_directory PROCESSED 
    File path where you would like the cleaned training and test sets to be stored.
  --train_expression_file_name TRAIN_EXPRESSION_FILE_NAME=${FILENAME_LEAD}_${SEED}_train_expression.RDS 
    Name for the training expression set output file
  --test_expression_file_name TEST_EXPRESSION_FILE_NAME=${FILENAME_LEAD}_${SEED}_test_expression.RDS
    Name for the test expression set output file
  --train_targets_file_name TRAIN_TARGETS_FILE_NAME=${FILENAME_LEAD}_${SEED}_train_targets.tsv
    Name for the training label set output file
  --test_targets_file_name TEST_TARGETS_FILE_NAME=${FILENAME_LEAD}_${SEED}_test_targets.tsv
    Name for the test label set output file
  --full_targets_file_name FULL_TARGETS_FILE_NAME=${FILENAME_LEAD}_${SEED}_full_targets.tsv
    Name for the combined training and test label set output file
  --seed $SEED
    seed for random number generation
  --train_percent $TRAIN_PERCENT 
    controls training/test partition

```

### 02-train_elasticnet

This script builds an elastic net logistic regression models from the training set and known target values set
that were cleaned and saved in script 01-clean_split_data.R.  The models are generated by a loop over the values
of the TRANSCRIPT_TAIL_ARRAY argument in the USER_ARGUMENTS section of run-sex-prediction-from-RNASeq.sh.


**Argument descriptions**
```
Values for all arguments are computed in run-sex-prediction-from-RNASeq.sh using values from the 
USER-SPECIFIED-ARGUMENTS and arguments computed for 01-clean_split_data.R.  i is the index of the
loop over the values of TRANSCRIPT_TAIL_ARRAY.

 --train_expression_file_name ${PROCESSED}/$TRAIN_EXPRESSION_FILE_NAME 
 --train_targets_file_name ${PROCESSED}/$TRAIN_TARGETS_FILE_NAME 
 --output_directory $MODELS 
 --model_object_file_name $MODEL_OBJECT_FILE_NAME 
 --model_transcripts_file_name $MODEL_TRANSCRIPTS_FILE_NAME 
 --model_coefs_file_name $MODEL_COEFS_FILE_NAME 
 --train_target_column $TRAIN_TARGET_COLUMN 
 --transcript_tail_percent ${TRANSCRIPT_TAIL_PERCENT_ARRAY[i]}
 
```
### 03-evaluate_model

This script evaluates the elastic net logistic regression models that were trained in script 
02-train_elasticnet.R from a training set and known target values set that was cleaned and saved 
in script 01-clean_split_data.R.  This script runs within the same loop described for script
02-train_elasticnet.R, and within a second nested loop over the values of the targetColumns
argument in the USER_SPECIFIED_ARGUMENTS section of run-sex-prediction-from-RNASeq.sh.  Output
files for each value of targetColumns go to a separate directory named with the current
targetColumns value.

Model evaluations are done with the caret::confusionMatrix and caret::twoClassSummary functions.
twoClassSummary can fail if one of the Confusion Matrix rows is all zeros.  In that case, no
twoClassSummary file is output for the transcript_tail_percent value in question.


**Argument descriptions**
```
Values for all arguments are computed in run-sex-prediction-from-RNASeq.sh using values from the 
USER-SPECIFIED-ARGUMENTS and arguments computed for 01-clean_split_data.R and 02-train_elasticnet.R.
i is the index of the loop over the values of TRANSCRIPT_TAIL_ARRAY; t is the index of the loop over
the values of targetColumns.

      --test_expression_file_name ${PROCESSED}/$TEST_EXPRESSION_FILE_NAME 
      --test_targets_file_name ${PROCESSED}/$TEST_TARGETS_FILE_NAME 
      --model_object_file_name ${MODELS}/$MODEL_OBJECT_FILE_NAME 
      --model_transcripts_file_name ${MODELS}/$MODEL_TRANSCRIPTS_FILE_NAME 
      --test_target_column ${t} 
      --output_directory $RESULTS_OUTPUT_DIRECTORY=${RESULTS}/${TEST_TARGET_COLUMN}
      --cm_set_file_name $CM_SET_FILE=${FILENAME_LEAD}_${SEED}_${TRANSCRIPT_TAIL_PERCENT_ARRAY[i]}_prediction_details.tsv
        A file containing prediction proabilities for each sample in the test set.
      --cm_file_name $CM_SET=${FILENAME_LEAD}_${SEED}_${TRANSCRIPT_TAIL_PERCENT_ARRAY[i]}_confusion_matrix.RDS
        A file containing the confusion matrix object out of caret::confusionMatrix.
      --summary_file_name $SUMMARY_FILE=${FILENAME_LEAD}_${SEED}_${TRANSCRIPT_TAIL_PERCENT_ARRAY[i]}_two_class_summary.RDS
        A file containing the two class summary object out of caret::twoClassSummary.

```

### 04-present_results

A Notebook that produces tables and plots from the output files of 02-train_elasticnet.R and 03-evaluate_model.R.  
Plots produced are Strength of Calls at Maximum Accuracy, Predictive Accuracy vs. Number of Training Transcripts,
Number of Non-Zero Features vs. Number of Training Transcripts.  Tables are Confusion Matrix at Maximum Accuracy
and Two Class Summary at Maximum Accuracy. Each of the tables and charts is produced for each value of the
targetColumns array specified in the USER-SPECIFIED ARGUMENTS section of run-sex-prediction-from-RNASeq.sh.


**Argument descriptions**
```
Values for all arguments are computed in run-sex-prediction-from-RNASeq.sh using values from the 
USER-SPECIFIED-ARGUMENTS and arguments computed for 01-clean_split_data.R, 02-train_elasticnet.R
and 03-evaluate_model.R.

  --results_dir=${RESULTS} 
  --model_dir=$MODELS 
  --cm_set=$CM_SET 
  --seed=$SEED 
  --target_columns=$targetColumns_to_pass
    A string constructed from the values of the targetColumns array in USER-SPECIFIED ARGUMENTS.

```
