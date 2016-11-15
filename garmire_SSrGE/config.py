"""
CONFIG file for SSrGE

Principal default values for SSrGE class parameters

The config file gives also parameters in order to
extract SNV and GE matrices from a given project

"""

######## SSrGE VARIABLE ##############################
TIME_LIMIT = 5
# time limit for one linear regression model
MIN_OBS_FOR_REGRESS = 10
# Min number of cell having non null gene expression
# to infer a sparse linear model
NB_THREADS = 4
# Number of threads to run in parallel
CROSS_VAL_NFOLD = 5
# Number of folds to perform the cross validation
######################################################
