

/*All methods discussed in SUW (2023)-Kappa paper
#Use gmm for the standard errors
#methods: tau_t,norm, tau,a_10
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
	
	
			tempvar troot 
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

	tempvar ips ipsxb numhat kappaw kappa_0 kappa_1 num1hat num0hat 
	tempname bips nums kappa_1s kappa_0s kappas late_a late_a1 late_a0 num1hats num0hats late_a10
	
quietly {
*Estimation of IPS

*If Logit:
	if "`zmodel'" == "logit" {
  `zmodel' `zvar' `xvarsips' if `touse'==1 
  matrix `bips' = e(b)
  predict double `ips'
 }
 else if "`zmodel'" == "probit" {

  `zmodel' `zvar' `xvarsips' if `touse'==1
  matrix `bips' = e(b)
  predict double `ips'
  
}
*If CBPS 
else if "`zmodel'" == "cbps" {

  `zmodel' `zvar' `xvarsips' if `touse'==1
 
    matrix `bips' = e(b)
    predict double `ipsxb'
    ge `ips' = logistic(`ipsxb')
  
}

 gen `numhat' = ((`zvar')/`ips')*`yvar'-(((1-`zvar'))/(1-`ips'))*`yvar' if `touse'==1


 gen `kappa_1' = ((`zvar')/`ips')*`tvar'-(((1-`zvar'))/(1-`ips'))*`tvar' if `touse'==1

 gen `kappa_0' = (1-`tvar')*((1-`zvar')-(1-`ips'))/(`ips'*(1-`ips')) if `touse'==1

 gen `kappaw' = 1- (`tvar'*(1-`zvar'))/(1-`ips')-((1-`tvar')*`zvar')/`ips' if `touse'==1


ge `num1hat' = `kappa_1'*`yvar' if `touse'==1

ge `num0hat' = `kappa_0'*`yvar' if `touse'==1


 sum `numhat' if `touse'==1 
matrix `nums' = r(mean)

 sum `kappa_1'  if `touse'==1 
matrix `kappa_1s' = r(mean)

 sum `kappa_0' if `touse'==1  
matrix `kappa_0s' = r(mean)

 sum `kappaw' if `touse'==1  
matrix `kappas' = r(mean)

 sum `num1hat' if `touse'==1  
matrix `num1hats' = r(mean)

 sum `num0hat' if `touse'==1  
matrix `num0hats' = r(mean)

matrix `late_a' =`nums'*invsym(`kappas')
matrix `late_a1' =`nums'*invsym(`kappa_1s')
matrix `late_a0' =`nums'*invsym(`kappa_0s')
matrix `late_a10' =`num1hats'*invsym(`kappa_1s')-`num0hats'*invsym(`kappa_0s')

matrix list `late_a'
matrix list `late_a1'
matrix list `late_a0'
matrix list `late_a10'

tempvar y1hat y0hat d1hat d0hat 
tempname by1 by0 bd1 bd0 denom1s denom0s num1s num0s num_norms denom_norms late_norm
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
matrix `late_norm' =`num_norms'*invsym(`denom_norms')
	
	matrix list `late_norm'


	
	//Moment Conditions:
	
	//1- IPS
	if "`zmodel'" == "logit" {
local eqips (eqips:`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))

local eqips_inst instruments(eqips: `xvarsips')

}
else if "`zmodel'" == "probit" {
local eqips (eqips:   ( (`zvar' - normal({zhat: `xvarsips' _cons}))/(normal({zhat:})*(1-normal({zhat:}))))*normalden({zhat:}))
local eqips_inst instruments(eqips: `xvarsips')

}
else if "`zmodel'" == "cbps" {
local eqips (eqips: (`zvar' - exp({zhat: `xvarsips' _cons})/(1+exp({zhat:})))/((exp({zhat:})/(1+exp({zhat:})))*(1/(1+exp({zhat:})))))

local eqips_inst instruments(eqips: `xvarsips')
}

//Other Moments
//psi_Delta
local eq_delta (eq_delta: ((`zvar'*`yvar')/(exp({zhat:})/(1+exp({zhat:})))-((1-`zvar')*`yvar')/(1-(exp({zhat:})/(1+exp({zhat:}))))- {deltap}))    

local eq_delta_inst instruments(eq_delta: )


//psi_Gamma
local eq_gamma (eq_gamma: (1-((1-`zvar')*`tvar')/(1-(exp({zhat:})/(1+exp({zhat:}))))- (`zvar'*(1-`tvar'))/(exp({zhat:})/(1+exp({zhat:})))- {gammap}))    

local eq_gamma_inst instruments(eq_gamma: )


//Psi_Gamma1
local eq_gamma1 (eq_gamma1: ((`zvar'*`tvar')/(exp({zhat:})/(1+exp({zhat:})))-((1-`zvar')*`tvar')/(1-(exp({zhat:})/(1+exp({zhat:}))))-{gamma1}))
local eq_gamma1_inst instruments(eq_gamma1: )


//psi_Gamma0
local eq_gamma0 (eq_gamma0: ((`zvar'*(`tvar'-1))/(exp({zhat:})/(1+exp({zhat:})))-((1-`zvar')*(`tvar'-1))/(1-(exp({zhat:})/(1+exp({zhat:}))))-{gamma0}))
local eq_gamma0_inst instruments(eq_gamma0: )


//Psi_Delta1
local eq_delta1 (eq_delta1:  (`tvar'*((`zvar'-(exp({zhat:})/(1+exp({zhat:}))))/((exp({zhat:})/(1+exp({zhat:})))*(1-(exp({zhat:})/(1+exp({zhat:}))))))*`yvar'-{delta1}))
local eq_delta1_inst instruments(eq_delta1: )

//Psi_Delta0
local eq_delta0 (eq_delta0:  ((1-`tvar')*(((1-`zvar')-(1-(exp({zhat:})/(1+exp({zhat:})))))/((exp({zhat:})/(1+exp({zhat:})))*(1-(exp({zhat:})/(1+exp({zhat:}))))))*`yvar'-{delta0}))
local eq_delta0_inst instruments(eq_delta0: )

//psi_mu1
local eq_mu1 (eq_mu1: ((`zvar'*(`yvar'-{mu1}))/(exp({zhat:})/(1+exp({zhat:}))))) 
local eq_mu1_inst instruments(eq_mu1: )

//psi_mu0
local eq_mu0 (eq_mu0: (((1-`zvar')*(`yvar'-{mu0}))/(1-(exp({zhat:})/(1+exp({zhat:})))))) 
local eq_mu0_inst instruments(eq_mu0: )

//psi_m1
local eq_m1 (eq_m1: ((`zvar'*(`tvar'-{m1}))/(exp({zhat:})/(1+exp({zhat:}))))) 
local eq_m1_inst instruments(eq_m1: )
//psi_m0
local eq_m0 (eq_m0: (((1-`zvar')*(`tvar'-{m0}))/(1-(exp({zhat:})/(1+exp({zhat:})))))) 
local eq_m0_inst instruments(eq_m0: )

//psi_taua
local eq_tau_a (eq_tau_a: ({tau_a} - {deltap}/{gammap}))

//psi_taua1
local eq_tau_a1 (eq_tau_a1: ({tau_a1} - {deltap}/{gamma1}))

//psi_taua0
local eq_tau_a0 (eq_tau_a0: ({tau_a0} - {deltap}/{gamma0}))

//psi_taua10
local eq_tau_a10 (eq_tau_a10: ({tau_a10} - ({delta1}/{gamma1}-{delta0}/{gamma0})))

//psi_taunorm

if `dmeanz0' !=0 & `dmeanz1' !=1  {
local eq_tau_norm (eq_tau_norm: ({tau_norm} - (({mu1}-{mu0})/({m1}-{m0}))))
}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
local eq_tau_norm (eq_tau_norm: ({tau_norm} - (({mu1}-{mu0})/({m1}))))	
}
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
	local eq_tau_norm (eq_tau_norm: ({tau_norm} - (({mu1}-{mu0})/(1-{m0}))))
}
tempname tau_a vc_tau_a r_tau_a var_tau_a tau_a1 vc_tau_a1 r_tau_a1 var_tau_a1 tau_a0 vc_tau_a0 r_tau_a0 var_tau_a0
tempname tau_norm vc_tau_norm r_tau_norm var_tau_norm tau_a10 vc_tau_a10 r_tau_a10 var_tau_a10
//Estimation:
//tau_a
matrix initial = (`bips', `nums', `kappas', `late_a')	
local k = colsof(initial )
matrix I = I(`k')  
gmm   `eqips' `eq_delta' `eq_gamma'  `eq_tau_a' if `touse',  `eqips_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce') iterate(0)
scalar `tau_a' = _b[tau_a:_cons]
matrix `vc_tau_a'= e(V)
scalar `r_tau_a' = rownumb(`vc_tau_a', "tau_a:_cons")
scalar `var_tau_a' = `vc_tau_a'[`r_tau_a', `r_tau_a'] 



//tau_a1
matrix initial = (`bips', `nums', `kappa_1s', `late_a1')	
local k = colsof(initial )
matrix I = I(`k')  
gmm   `eqips' `eq_delta' `eq_gamma1'  `eq_tau_a1' if `touse',  `eqips_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce') iterate(0)
scalar `tau_a1' = _b[tau_a1:_cons]
matrix `vc_tau_a1'= e(V)
scalar `r_tau_a1' = rownumb(`vc_tau_a1', "tau_a1:_cons")
scalar `var_tau_a1' = `vc_tau_a1'[`r_tau_a1', `r_tau_a1'] 



//tau_a0
matrix initial = (`bips', `nums', `kappa_0s', `late_a0')	
local k = colsof(initial )
matrix I = I(`k')  
gmm   `eqips' `eq_delta' `eq_gamma0'  `eq_tau_a0' if `touse',  `eqips_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce') iterate(0)
scalar `tau_a0' = _b[tau_a0:_cons]
matrix `vc_tau_a0'= e(V)
scalar `r_tau_a0' = rownumb(`vc_tau_a0', "tau_a0:_cons")
scalar `var_tau_a0' = `vc_tau_a0'[`r_tau_a0', `r_tau_a0'] 


//tau_a10
matrix initial = (`bips', `num1hats', `kappa_1s', `num0hats', `kappa_0s', `late_a10')	
local k = colsof(initial )
matrix I = I(`k')  
gmm   `eqips' `eq_delta1' `eq_gamma1' `eq_delta0' `eq_gamma0'  `eq_tau_a10' if `touse',  `eqips_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce') iterate(0)
scalar `tau_a10' = _b[tau_a10:_cons]
matrix `vc_tau_a10'= e(V)
scalar `r_tau_a10' = rownumb(`vc_tau_a10', "tau_a10:_cons")
scalar `var_tau_a10' = `vc_tau_a10'[`r_tau_a10', `r_tau_a10'] 


//tau_norm

if `dmeanz0' !=0 & `dmeanz1' !=1  {
matrix initial = (`bips', `num1s', `num0s', `denom1s', `denom0s', `late_norm')	
local k = colsof(initial )
matrix I = I(`k')  
gmm  `eqips' `eq_mu1' `eq_mu0'  `eq_m1' `eq_m0'  `eq_tau_norm' if `touse' ,  `eqips_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce') iterate(0)

}
else if `dmeanz0' == 0 & `dmeanz1' !=1 {
matrix initial = (`bips', `num1s', `num0s', `denom1s', `late_norm')	
local k = colsof(initial )
matrix I = I(`k')  
gmm  `eqips' `eq_mu1' `eq_mu0'  `eq_m1' `eq_tau_norm' if `touse' ,  `eqips_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce') iterate(0)

}
else if `dmeanz0' !=0 & `dmeanz1' ==1 {
matrix initial = (`bips', `num1s', `num0s', `denom0s', `late_norm')	
local k = colsof(initial )
matrix I = I(`k')  
gmm  `eqips' `eq_mu1' `eq_mu0' `eq_m0'  `eq_tau_norm' if `touse' ,  `eqips_inst' onestep winitial(I) from(initial) quickderivatives vce(`vce') iterate(0)
	
}
scalar `tau_norm' = _b[tau_norm:_cons]
matrix `vc_tau_norm'= e(V)
scalar `r_tau_norm' = rownumb(`vc_tau_norm', "tau_norm:_cons")
scalar `var_tau_norm' = `vc_tau_norm'[`r_tau_norm', `r_tau_norm'] 
 }

	tempname b V
	matrix `b' = (`tau_a', `tau_a1', `tau_a0', `tau_a10', `tau_norm')
	matrix `V' = (`var_tau_a', 0, 0, 0, 0\ 0, `var_tau_a1', 0,0,0 \ 0, 0,`var_tau_a0', 0,0\ 0,0,0, `var_tau_a10', 0\  0, 0, 0, 0, `var_tau_norm' )
	matrix rownames `b' = " "
	matrix colnames `b' = "tau_a" "tau_a1" "tau_a0" "tau_a10" "tau_norm"
	matrix rownames `V' = "tau_a" "tau_a1" "tau_a0" "tau_a10" "tau_norm"
	matrix colnames `V' = "tau_a" "tau_a1" "tau_a0" "tau_a10" "tau_norm"
 
 
 	di as text "LATE Estimation"
	di as text "Outcome        : "  "`yvar'"
	di as text "Treatment      : "  "`tvar'"
	di as text "Instrument     : "  "`zvar'"
	if "`zmodel'" == "logit" {
	di as text "IPS Model      : logit"
	}
		else if "`zmodel'" == "probit" {
	di as text "IPS Model      : probit"
	}
	else if "`zmodel'" == "cbps" {
	di as text "IPS Model      : CBPS"
	}
	di

	ereturn post `b' `V', esample(`touse')
	ereturn local cmd "kappalate"
	ereturn local depvar `yvar'
	ereturn local tvar `tvar'
	ereturn local tvar `zvar'
	ereturn local stat `stat'
	ereturn local zmodel `zmodel'
	ereturn local vcetype "Robust"
	ereturn local vce `vce'
	ereturn local cmdline `"late `0'"'
	ereturn local title "LATE Estimation"
    


	ereturn display
	
end
	
	
