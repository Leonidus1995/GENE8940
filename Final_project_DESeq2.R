library(tximportData)

# Variable dir stores the complete path to the directory
#extdata within the tximportData package
dir <- system.file("extdata", package = "tximportData")
list.files(dir)

# Next, we create a named vector pointing to the quantification files
