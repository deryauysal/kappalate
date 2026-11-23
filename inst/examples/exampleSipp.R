# use the SIPP data
# install.packages("haven")  # if needed
library(kappalate)
library(haven)

# SIPP data shipped with the package
data <- as.data.frame(read_dta(system.file("extdata", "sipp.dta", package = "kappalate")))


# Drop rows where kwage is missing (NA) or educ is missing (NA) or rsncode == 999
data <- data[!(is.na(data$kwage) | is.na(data$educ) | data$rsncode == 999), ]

# Create a new variable lwage as the natural log of kwage
data$lwage <- log(data$kwage)


kappalate_cbps <- kappalate(lwage ~ age_5 | nvstat | rsncode  , data = data, zmodel = "cbps", std = "on",  which = "all")
kappalate_logit <- kappalate(lwage ~ age_5 | nvstat | rsncode  , data = data, zmodel = "logit", std = "on",  which = "all")

print(kappalate_cbps)
print(kappalate_logit)