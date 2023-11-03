        use https://economics.mit.edu/sites/default/files/inline-files/sipp2.dta, clear

        drop if kwage==. | educ==. | rsncode==999

        generate double lwage = ln(kwage)

        kappalate lwage (nvstat=rsncode) age_5, zmodel(cbps) which(norm)

        kappalate lwage (nvstat=rsncode) age_5, zmodel(cbps) which(all)

        kappalate lwage (nvstat=rsncode) age_5, zmodel(logit) which(all)
