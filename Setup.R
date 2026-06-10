## Packages
library(tidyverse)
## Data Input
# Read the data from the rdata.file and summarize the data
load(file.choose())
data <- data.frame(MultinomialExample)
summary(data)
# There are 500 observations in the dataset, so we randomly choose 400 observations for the training and 100 for the testing 
## We randomly separate the data into 3 groups then test the factor analyzer
set.seed(10000)
train_idx <- sample(1:nrow(data), size = 400, replace = FALSE)
training <- data[train_idx, ]
testing  <- data[-train_idx, ]
group <- sample(1:3, size = nrow(training), replace = TRUE)
training$group <- group
for (i in 1:3) {
  assign(paste0("group", i), 
         training %>% filter(group == i) %>% select(-group))
  cat("Group", i, ":", nrow(get(paste0("group", i))), "observations\n")
}
covariate_matrix <- function(df) {
  x_cols <- grep("^x", names(df), value = TRUE)
  return(as.matrix(df[, x_cols]))
}
for (i in 1:3) {
  assign(paste0("X", i), 
         covariate_matrix(get(paste0("group", i))))
}
