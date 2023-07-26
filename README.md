# QHEBO-paper
Source files and code used for the paper "Health outcomes in Bulgaria: simulated effects of obesogenic environmental changes in adulthood versus childhood"


The "0 Create BMI matrices" Rmd file needs to be run first. It does two things:
  1) Extrapolates the BMI categories for men and women seperately beyond 2016. Saves them in the main directory as "BMIdfs.RData".
  2) Creates transition matrices between the three BMI states. Saves them in the main directory as "P1matrix_SEX_YEAR.RData".

Then "1 Main.RData" can work properly.
