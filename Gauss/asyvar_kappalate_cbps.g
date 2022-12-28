//asymptotic variance for normalized tau


proc(1)= asyvar_kappalate(struct PV p, struct DS ds_struct,k); 
local indep, dep, treat, inst, N, pf ;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   inst = ds_struct[4].dataMatrix;  
   pf = ds_struct[5].Type; //probit(2) or logit(1) for the IPS 
   N = ds_struct[6].Type;
 local alp, gammah, deltah, tauh;
   
     alp = pvUnpack(p, "alphahat");
     gammah = pvUnpack(p, "gammahat");
     deltah = pvUnpack(p, "deltahat");
     tauh = pvUnpack(p, "tauhat");
    
local psest, dpsest, Edeltaalp, Egammaalp, EH, psialp, psidelta, psigamma,delpsi;


if pf==1;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psialp = (indep.*(inst-psest)./(psest.*(1-psest))); //Prob. function logit.
delpsi = psest.*(1-psest).*indep;

endif;




Edeltaalp = - meanc((((dep.*inst)./(psest.^2))+((dep.*(1-inst))./((1-psest).^2))).*delpsi);

if k==0; //kappa0
//Egammaalp = - meanc(((((1-treat).*inst)./(psest.^2))+(((1-treat).*(1-inst))./((1-psest).^2))).*psest.*(1-psest).*indep);
Egammaalp= -meanc(((((treat-1).*inst)./(psest.^2))+(((treat-1).*(1-inst))./((1-psest).^2))).*delpsi);    
elseif k==1; // kappa1
Egammaalp = - meanc((((treat.*inst)./(psest.^2))+((treat.*(1-inst))./((1-psest).^2))).*delpsi);

elseif k==2; // kappa    
//Egammaalp = - meanc(((((1-treat).*inst)./(psest.^2))+(((treat).*(1-inst))./((1-psest).^2))).*psest.*(1-psest).*indep);
Egammaalp =  meanc(((((1-treat).*inst)./(psest.^2))-((treat.*(1-inst))./((1-psest).^2))).*delpsi);
   
endif;

EH = -(1/N)*psialp'*psialp;

psidelta= ((inst./psest).*dep-((1-inst)./(1-psest)).*dep)-deltah;
if  k==0;
psigamma = ((inst./psest).*(treat-1)-((1-inst)./(1-psest)).*(treat-1))-gammah;
//psigamma = (((1-inst)./(1-psest))-(inst./psest)+(inst./psest).*treat-((1-inst)./(1-psest)).*treat)-gammah;
elseif k==1;
  psigamma = ((inst./psest).*treat-((1-inst)./(1-psest)).*treat)-gammah;
elseif k==2;
psigamma = 1-(treat.*(1-inst))./(1-psest)-((1-treat).*inst)./psest-gammah;
endif;



local nix1, nix2, Vnix;
nix1=-((1/gammah).*(Edeltaalp)-(deltah/gammah^2).*(Egammaalp))'*inv(-EH)*((1/gammah).*(Edeltaalp)-(deltah/gammah^2).*(Egammaalp));
nix2=meanc(((1/gammah).*(psidelta)-(deltah/gammah^2).*(psigamma)).^2);
Vnix=nix1+nix2;

retp(sqrt((1/N)*Vnix));
endp;

