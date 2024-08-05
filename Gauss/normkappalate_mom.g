
@
Moment functions for tau_a,10
@


proc normkappalate_mom(struct PV p, struct DS ds_struct); 
local indep, dep, treat, inst, pf, N ;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   inst = ds_struct[4].dataMatrix;
     pf = ds_struct[5].Type; //probit or logit for the IPS 
      N = ds_struct[6].Type;
local alp, psi1,psi2,psi3, psi4,psi5,psi6, psest,dpsest, delta1h, delta0h, gamma1h, gamma0h, tauh;
     alp = pvUnpack(p, "alphahat");
     delta1h = pvUnpack(p, "delta1hat");
     gamma1h = pvUnpack(p, "gamma1hat");
     delta0h = pvUnpack(p, "delta0hat");
     gamma0h = pvUnpack(p, "gamma0hat");
     tauh = pvUnpack(p, "tauhat");


if pf==1;//if logit
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = meanc(indep.*(inst-psest)); //Prob. function logit.
elseif pf==2;//if probit
psest=cdfn((indep)*alp); //Prob. function probit.
dpsest=pdfn((indep)*alp); //density. function probit.
psi1 = meanc(indep.*dpsest.*(inst-psest)./(psest.*(1-psest))); 
elseif pf==3;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = meanc(indep.*(inst-psest)./(psest.*(1-psest))); //Prob. function logit.
endif;


psi2 = meanc(treat.*((inst-psest)./(psest.*(1-psest))).*dep-delta1h);
psi3 = meanc(treat.*((inst-psest)./(psest.*(1-psest))) - gamma1h);
psi4 = meanc((1-treat).*(((1-inst)-(1-psest))./(psest.*(1-psest))).*dep- delta0h);
psi5 = meanc((1-treat).*(((1-inst)-(1-psest))./(psest.*(1-psest)))- gamma0h);
psi6 = (delta1h/gamma1h-delta0h/gamma0h-tauh);    


retp (psi1|psi2|psi3|psi4|psi5|psi6);

endp;

