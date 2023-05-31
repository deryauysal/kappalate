clear

do "kappalate.ado"
do "cbps.ado"
cls 
clear matrix
clear mata
set maxvar 120000

use "Abadie2003_data.dta"

egen incst= std(inc)
egen age25st= std(age25)
egen age25sqst= std(age25sq)
egen fsizest = std(fsize)

kappalate net_tfa p401 e401 incst age25st age25sqst marr fsizest

kappalate net_tfa p401 e401 incst age25st age25sqst marr fsizest, zmodel(cbps)

kappalate net_tfa p401 e401 incst age25st age25sqst marr fsizest, zmodel(probit)
