


//Collect all moments for weighting here as functions of:
//D: Dependent variable 
//Z: Treatment
//X: Covariate matrix
//Ps specification-> logit

proc normlate_mom_cbps(struct PV p, struct DS ds_struct); 
local indep, dep, treat, inst, pf, N ;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   inst = ds_struct[4].dataMatrix;    
       pf = ds_struct[5].Type; //probit or logit for the IPS 
      N = ds_struct[6].Type;
local alp, psi1,psi2,psi3, psi4,psi5,psi6,psi7, psi8, psest, dpsest, mu1, mu0, mu1d, mu0d, gammah, deltah, tauh;

     alp = pvUnpack(p, "alphahat");
     mu1 = pvUnpack(p, "mu1hat");
     mu0 = pvUnpack(p, "mu0hat");
     mu1d = pvUnpack(p, "m1hat");
     mu0d = pvUnpack(p, "m0hat");
     deltah = pvUnpack(p, "deltahat");
     gammah = pvUnpack(p, "gammahat");
     tauh = pvUnpack(p, "tauhat");

if pf==1;//if logit
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = meanc(indep.*(inst-psest)./(psest.*(1-psest))); //Prob. function logit.
endif;


psi2=meanc(inst.*(dep-mu1)./psest);
psi3=meanc((1-inst).*(dep-mu0)./(1-psest));
psi4=meanc(inst.*(treat-mu1d)./psest);
psi5=meanc((1-inst).*(treat-mu0d)./(1-psest));
psi6=mu1-mu0-deltah;
psi7=mu1d-mu0d-gammah;
psi8=deltah/gammah-tauh;

retp (psi1|psi2|psi3|psi4|psi5|psi6|psi7|psi8);

endp;
