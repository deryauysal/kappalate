


//Collect all moments for weighting here as functions of:
//D: Dependent variable 
//Z: Treatment
//X: Covariate matrix
//Ps specification-> logit

proc kappalate_mom(struct PV p, struct DS ds_struct); 
local indep, dep, treat, inst, pf, N ;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   inst = ds_struct[4].dataMatrix;
   pf = ds_struct[5].Type; //probit or logit for the IPS 
   N = ds_struct[6].Type;    
local alp, psi1,psi2,psi3, psi4, tauh, psest, dpsest, gammah, deltah;
     alp = pvUnpack(p, "alphahat");
     gammah = pvUnpack(p, "gammahat");
     deltah = pvUnpack(p, "deltahat");
     tauh = pvUnpack(p, "tauhat");

if pf==1;//if logit
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = meanc(indep.*(inst-psest)./(psest.*(1-psest))); //Prob. function logit.
endif;


psi2=meanc((inst./psest).*dep-((1-inst)./(1-psest)).*dep)-deltah;
//psi3=meanc(1-(inst./psest)+(inst./psest).*treat-((1-inst)./(1-psest)).*treat)-gammah;
psi3=meanc(1-(treat.*(1-inst))./(1-psest)-((1-treat).*inst)./psest-gammah);
psi4=deltah/gammah-tauh;

retp (psi1|psi2|psi3|psi4);

endp;
