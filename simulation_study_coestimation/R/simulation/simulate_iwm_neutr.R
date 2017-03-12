#!/usr/bin/Rscript --no-site-file --no-save

library(coala)
library(jaatha)
library(testJaatha)

# Load the model
model_iwm_neutr <- NULL
load("cache/models.Rda")
stopifnot(!is.null(model_iwm_neutr))

# Load the test dataset
test_data_file <- "data/model_iwm.Rda"
if (!file.exists(test_data_file)) stop("Test data not found")
load(test_data_file)

# Execute the jaatha analyses
testJaatha(dm = model_iwm_neutr, test_data = test_data,
           folder = "data/iwm_neutr", seed = 2202,
           sim = 160, init_method = "zoom-in")
