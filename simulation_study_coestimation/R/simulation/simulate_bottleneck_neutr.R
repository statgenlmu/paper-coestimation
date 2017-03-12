#!/usr/bin/Rscript --no-site-file --no-save

library(coala)
library(jaatha)
library(testJaatha)

# Load the model
model_bn_neutr <- NULL
load("cache/models.Rda")
stopifnot(!is.null(model_bn_neutr))

# Load the test dataset
test_data_file <- "data/model_bottleneck.Rda"
if (!file.exists(test_data_file)) stop("Test data not found")
load(test_data_file)

# Execute the jaatha analyses
testJaatha(dm = model_bn_neutr, test_data = test_data,
           folder = "data/bottleneck_neutr", seed = 6606,
           sim = 160, init_method = "zoom-in")
