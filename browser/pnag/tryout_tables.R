library(dplyr)
library(kableExtra)

# make dummy table wit 8 cols and 10 rows
df <- data.frame(matrix(rnorm(80), nrow = 10))

# use kableExtra to make table
