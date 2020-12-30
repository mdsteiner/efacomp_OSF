##### Compare original results from SPSS 23 and 27

if (!require(EFAtools)) install.packages("EFAtools"); library(EFAtools)

# Print decimal instead of exponential notation of numbers
options(scipen = 999)

# # prepare data
SPSS_23$NEO <- NULL

#### Compare using the Compare function

paf_SPSS_23_27 <- list()
var_SPSS_23_27 <- list()
pro_SPSS_23_27 <- list()

for(i in seq_along(SPSS_23)){
  paf_SPSS_23_27[[i]] <- COMPARE(SPSS_23[[i]]$paf_load, SPSS_27[[i]]$paf_load)
  var_SPSS_23_27[[i]] <- COMPARE(SPSS_23[[i]]$var_load, SPSS_27[[i]]$var_load)
  pro_SPSS_23_27[[i]] <- COMPARE(SPSS_23[[i]]$pro_load, SPSS_27[[i]]$pro_load)
}

names(paf_SPSS_23_27) <- names(SPSS_27)
names(var_SPSS_23_27) <- names(SPSS_27)
names(pro_SPSS_23_27) <- names(SPSS_27)