### ** Examples

# use the SIPP data
install.packages("haven")  # Install if you haven't
library(haven)

data <- read_dta("~/Documents/GitHub/kappalate/stata/sipp.dta")  # Replace with your actual file path


# Drop rows where kwage is missing (NA) or educ is missing (NA) or rsncode == 999
data <- data[!(is.na(data$kwage) | is.na(data$educ) | data$rsncode == 999), ]

# Create a new variable lwage as the natural log of kwage
data$lwage <- log(data$kwage)


kappalate_cbps <- kappalate(lwage ~ age_5 | nvstat | rsncode  , data = data, zmodel = "cbps", std = "on",  which = "all")
kappalate_logit <- kappalate(lwage ~ age_5 | nvstat | rsncode  , data = data, zmodel = "logit", std = "on",  which = "all")

print(kappalate_cbps)
print(kappalate_logit)