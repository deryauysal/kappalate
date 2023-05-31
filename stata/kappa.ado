

/*All methods discussed in SUW (2023)-Kappa paper
#Use gmm for the standard errors
#methods: tau_cb, tau_t,norm, tau,a_10
tau_a, tau_a,1 tau_a,0
*/


capture program drop kappalate

program kappalate, eclass 

version 17 

syntax  namelist [if] [,zmodel(string) vce(string)]
	
	local orig `0'
	tokenize `namelist'
	local yvar `1'  //outcome
	macro shift
	local tvar `1'  //treatment
	macro shift 
	local zvar `1'	 //instrument	
	macro shift
	local xvarsips `*'		//predetermined variables (if models use same set of covariates, otherwise empty)
	
	if "`zmodel'" == "" {
		local zmodel logit
	}
	
	if "`vce'" == "" {
		local vce robust
	}
	
	marksample touse
	
	
	/*		tempvar troot 
		tempname valuest
		qui tab `tvar' if `touse', gen(`troot') matrow(`valuest')
		if rowsof(`valuest') != 2 {
			noi di as error "The treatment variable must only take 2 values in the sample."
			error 450
			exit
		}
	    quietly summarize `tvar' if `touse'
	    if r(min)!=0 | r(max)!=1 {
		noi display as error "Treatment must only take on values zero or one"
		error 450
		exit
	    }
		
		tempvar zroot 
		tempname valuesz
		qui tab `zvar' if `touse', gen(`zroot') matrow(`valuesz')
		if rowsof(`valuesz') != 2 {
			noi di as error "The instrumental variable must only take 2 values in the sample."
			error 450
			exit
		}			
	    quietly summarize `zvar' if `touse'
	    if r(min)!=0 | r(max)!=1 {
		noi display as error "Instrument must only take on values zero or one"
		error 450
		exit
	    }	
		
		tempname dmeanz1 dmeanz0
//Check the compliance:
	    quietly summarize `tvar' if `zvar' == 1 & `touse'==1
		scalar `dmeanz1' = r(mean)
		quietly summarize `tvar' if `zvar'  == 0 & `touse'==1
		scalar `dmeanz0' = r(mean)
*/
	tempvar ips ipsxb numhat kappaw kappa_0 kappa_1 num1hat num0hat
	tempname bips nums kappa_1s kappa_0s kappas late_a late_a1 late_a0 
	

*Estimation of IPS

*If Logit:
	if "`zmodel'" == "logit" {
  `zmodel' `zvar' `xvarsips' if `touse'==1 [pw = `samplew']
  matrix `bips' = e(b)
  predict double `ips'
 }
 
*If CBPS 
else if "`zmodel'" == "cbps" {

  `zmodel' `zvar' `xvarsips' if `touse'==1, samplew(`samplew')
 
    matrix `bips' = e(b)
    predict double `ipsxb'
    ge `ips' = logistic(`ipsxb')
  
}
/*
qui gen `numhat' = ((`zvar')/`ips')*`yvar'-(((1-`zvar'))/(1-`ips'))*`yvar' if `touse'==1


qui gen `kappa_1' = ((`zvar')/`ips')*`tvar'-(((1-`zvar'))/(1-`ips'))*`tvar' if `touse'==1

qui gen `kappa_0' = (1-`tvar')*((1-`zvar')-(1-`ips'))/(`ips'*(1-`ips')) if `touse'==1

qui gen `kappaw' = 1- (`tvar'*(1-`zvar'))/(1-`ips')-((1-`tvar')*`zvar')/`ips' if `touse'==1


ge `num1hat' = `kappa_1'*`yvar' if `touse'==1

ge `num0hat' = `kappa_0'*`yvar' if `touse'==1


qui sum `numhat' if `touse'==1 
matrix `nums' = r(mean)

qui sum `kappa_1'  if `touse'==1 
matrix `kappa_1s' = r(mean)

qui sum `kappa_0' if `touse'==1  
matrix `kappa_0s' = r(mean)

qui sum `kappa' if `touse'==1  
matrix `kappas' = r(mean)

matrix `late_a' =`nums'*invsym(`kappas')
matrix `late_a1' =`nums'*invsym(`kappa_1s')
matrix `late_a0' =`nums'*invsym(`kappa_0s')
*/
/*
tempvar y1hat y0hat d1hat d0hat 
tempname by1 by0 bd1 bd0 denom1s denom0s num1s num0s num_norms denom_norms tau_norm
/*Normalized IPW-LATE*/
qui regress `yvar'  if `zvar'==1 & `touse'==1 [pw = 1/`ips']
qui matrix `by1' = e(b)
qui predict double `y1hat'

 
 
 qui regress `yvar'  if `zvar'==0 & `touse'==1 [pw = 1/(1-`ips')]
qui matrix `by0' = e(b)
qui predict double `y0hat'

//Conditional mean of d1 & d0#


	    if `dmeanz1'==1 {
        ge double `d1hat'= 1
	    }
		else {

qui regress `tvar'  if `zvar'==1  & `touse'==1  [pw = 1/`ips']
matrix `bd1' = e(b)
qui predict double `d1hat'
}
	    
	    if `dmeanz0'==0 {
        ge double `d0hat' = 0
	    }
		else {
qui regress `tvar'  if `zvar'==0  & `touse'==1 [pw =  1/(1-`ips')]
matrix `bd0' = e(b)
qui predict double `d0hat'
		}
		
		
qui sum `d1hat' if `touse'==1 
matrix `denom1s' = r(mean)
qui sum `d0hat'  if `touse'==1 
matrix `denom0s' = r(mean)
qui sum `y1hat' if `touse'==1  
matrix `num1s' = r(mean)
qui sum `y0hat' if `touse'==1 
matrix `num0s' = r(mean)
matrix `num_norms' = `num1s'-`num0s'
matrix `denom_norms' = `denom1s'-`denom0s'
matrix `tau_norm' =`num_norms'*invsym(`denom_norms')
	
	
	
	matrix list `tau_norm'*/
	
	/*matrix list `tau_a'
	matrix list `tau_a1'
	matrix list `tau_a0'
*/

	
end
	
	
