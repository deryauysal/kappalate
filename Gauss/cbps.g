@  Estiamtion of IPS (or PS) by CBPS based on Imai Ratkovic 2014

Derya Uysal

All rights reserved. Commercial use prohibited.

Using this programme for noncommercial purposes is free, if appropriate reference is given.

Use Gauss 22 

Version 1.0

        - cbps: procedure to get the cbps coeffecient estimates and se. uses as imput the moment function "cbps_mom"
        and  "cbps_crossmom"
        return:
   outcbps.pars=  estcoef; *coeffcients estimates of the probability model
   outcbps.parsse = se_coeff; *se based on M-estimation asymptotics
   outcbps.retcode = out.retcode; *return code to the non-linear equation solving
	    - cbps_mom: covariate balancing moment function
        - cbps_crossmom: sample counterpart of the middle term of the var-cov matrix
output:


	

@

@ ========================================================================================================= @
@ ======================== BEGIN of PROGRAMME ================================================================ @
@ ========================================================================================================= @

struct OutputsCBPS {
   matrix pars; // ips coeff, to means in numerator tow means in denominator
   matrix parsse;
   matrix retcode;
};



proc cbps(&f1,&f2, struct DS ds_struct);
    
    library dc;
local indep, dep, treat, inst, pf, N, v_mest, av_late, se_late;
   indep = ds_struct[1].dataMatrix; 
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   inst = ds_struct[4].dataMatrix; 
   pf = ds_struct[5].Type;
   N = ds_struct[6].Type;
     local f1:proc; //one of the 5 
    local f2:proc; // procedure for crossproduct of the moments   
        //starting values
    local alp_start;
    
      //alp_start=inv(indep'indep)*indep'inst;  
    // Step one: Declare dc control structure
struct dcControl dcCt;
 
// Initialize dc control structure
dcCt = dcControlCreate();
 
// Step two: Describe data names
dcSetYVar(&dcCt, inst);

// Name of independent variable
dcSetXVars(&dcCt, indep[.,2:cols(indep)]);

// Step three: Declare dcOut struct 
struct dcout dcout1;
 
// Step four: Call binary logit procedure
   
        dcout1 = binaryLogit(dcCt);
    alp_start=pvGetParVector(dcout1.par);
struct PV par;
par=pvCreate();
par= pvPack(par,alp_start, "alphahat");
local  estcoef,av_coef,se_coeff ;
struct eqSolvemtControl c;
c = eqSolvemtControlCreate();
c.maxIters=10000;    
c.output = 0;
// Declare output structure to hold results
struct eqSolvemtOut out;
// Solve the system of equations
out =  eqSolvemt(&f1, par, ds_struct, c);

    estcoef=pvGetParVector(out.par); //all coeffecient estimates

    
    v_mest=(1/N)*inv(gradMT(&f1,out.par, ds_struct))*f2(out.par, ds_struct)*inv(gradMT(&f1,out.par, ds_struct))';

    se_coeff=sqrt(diag(v_mest));
    //Define structure

   struct OutputsCBPS outcbps;
   outcbps.pars=  estcoef;
   outcbps.parsse = se_coeff;
   outcbps.retcode = out.retcode;
   retp(outcbps);
endp;
    



proc cbps_mom(struct PV p, struct DS ds_struct); 
local indep, dep, treat, inst, pf, N ;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   inst = ds_struct[4].dataMatrix; 
   pf = ds_struct[5].Type; //probit or logit for the IPS 
   N = ds_struct[6].Type;
local alp, psi1,psi2,psi3, psi4, psest, dpsest, gammah, deltah, tauh;
     alp = pvUnpack(p, "alphahat");


psest=exp((indep)*alp)./(1+exp((indep)*alp)); 
psi1 = meanc(indep.*(inst-psest)./(psest.*(1-psest))); 

retp (psi1);

endp;

 
proc cbps_crossmom(struct PV p, struct DS ds_struct); 
local indep, dep, treat, inst, pf, N ;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   inst = ds_struct[4].dataMatrix; 
   pf = ds_struct[5].Type; //probit or logit for the IPS 
   N = ds_struct[6].Type;
local alp, psi1,psi2,psi3, psi4, psest, dpsest, gammah, deltah, tauh;
     alp = pvUnpack(p, "alphahat");


psest=exp((indep)*alp)./(1+exp((indep)*alp)); 
psi1 = (indep.*(inst-psest)./(psest.*(1-psest))); 

retp ((1/N)*psi1'*psi1);

endp;   
