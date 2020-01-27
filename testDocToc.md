# Sex Prediction from RNASeq

The 01-clean_split_data.R, 02-train_elasticnet.R, 03-evaluate_model.R, 04-present_results.Rmd pipeline trains and evaluates an elasticnet logistic regression model to predict sex from RNASeq data.
The training features are gene expression transcripts, and the training labels are reported_gender values for each sample.

The pipeline is a response to issue [#84](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/84) which was raised to check whether, in some histologies, silencing might be breaking down, potentially resulting in changes in X inactivation.
Based on the accuracy achieved here, this is probably not happening, and this classifier can be helpful to predict values for datasets without annotated sex information.

## How to run sex-prediction-from-seq

To run the full pipeline of data prep, model building, model evaluation and results presentation, call the bash script:

```
bash run-sex-prediction-from-RNASeq.sh
```
This bash script requires arguments that are set in the USER-SPECIFIED ARGUMENTS section of the script.  
The script will return:

- Plots and tables in a notebook: `04-present_results.html`.  
  - Tables: confusion matrix and two class summary reports for the model with maximum accuracy.  
  - Plots: strength of calls on the test set for the maximum accuracy model, predictive accuracy vs. number of training transcripts and number non-zero features vs. number of training transcripts.
- Files organized in subfolders of the folder containing run-sex-prediction-from-RNASeq.sh as follows:
  - `processed_data` holds the output of script 01-clean_split_data.R.
  Five files: expression and target files for training and testing plus a file containing the training and testing target values combined.
  - `models` holds the output of script 02-train_elasticnet.R.
  Three files for each value in the run-sex-prediction-from-RNASeq.sh argument TRANSCRIPT_TAIL_PERCENT_ARRAY:
    - the best fitting model object for the given TRANSCRIPT_TAIL_PERCENT_ARRAY value, 
    - the indices of the training transcripts used in that model
    - the non-zero coefficients of the model.
  - `results` holds the output of script 03-evaluate_model.R.  
  Three files for each value in the run-sex-prediction-from-RNASeq.sh argument TRANSCRIPT_TAIL_PERCENT_ARRAY: 
    - the caret::confusionMatrix for the model for the given TRANSCRIPT_TAIL_PERCENT_ARRAY value, 
    - the caret::twoClassSummary object for the given TRANSCRIPT_TAIL_PERCENT_ARRAY value (Note: caret::twoClassSummary function can fail under certain conditions, so twoClassSummary files may be missing for certain TRANSCRIPT_TAIL_PERCENT_ARRAY values.)
    - prediction probabilities for each sample in the test set.
  - `results/plots` holds .png copies of the plots presented in the `04-present_results.html` notebook.
  
## Summary of Methods


Input data is split into training and test partitions according to the user-specified argument TRAIN_PERCENT.


The glmnet package is used to fit an elastic net logistic regression model via penalized maximum likelihood.
The glmnetUtils package is used to do elastic net cross-validation for alpha and lambda simultaneously.


The caret package is used to generate a confusion matrix object and a two-class summary object for each model.
Definitions of the statistics presented are found in the documentation for the [confusion matrix function](http://topepo.github.io/caret/measuring-performance.html#measures-for-predicted-classes) and the [two-class summary function](http://topepo.github.io/caret/measuring-performance.html#measures-for-class-probabilities).


## General usage of scripts

**Overall notes about these scripts:**
- The scripts are sequential as noted by their number.
- All file path-related arguments assume the file path specified is relative to `OpenPBTA-analysis/analyses/sex-prediction-from-RNASeq`.
- The scripts create user-specified output directories that do not exist at run time.
- The scripts add files whose names do not match existing files, if any, in the output directories.  
Output files whose names match existing files overwrite the previous versions.

