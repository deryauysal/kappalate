//asymptotic variance for normalized tau


proc(1)= asyvar_normlate(struct PV p, struct DS ds_struct); 
local indep, dep, treat, inst, N, pf ;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   inst = ds_struct[4].dataMatrix;  
   pf = ds_struct[5].Type; //probit or logit for the IPS 
   N = ds_struct[6].Type;

 local alp, mu1, mu0, mu1d, mu0d, gammah, deltah, tauh;
   
      alp = pvUnpack(p, "alphahat");
     mu1 = pvUnpack(p, "mu1hat");
     mu0 = pvUnpack(p, "mu0hat");
     mu1d = pvUnpack(p, "m1hat");
     mu0d = pvUnpack(p, "m0hat");
     deltah = pvUnpack(p, "deltahat");
     gammah = pvUnpack(p, "gammahat");
     tauh = pvUnpack(p, "tauhat");  
    
  local psest,dpsest, Emu1alp, Emu0alp, Em1alp, Em0alp, psimu1,psimu0, psim1,psim0,psialp, EH, delpsi;
    if pf==1;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psialp = (indep.*(inst-psest)); //Prob. function logit.
delpsi = psest.*(1-psest).*indep;    
elseif pf==2;
psest=cdfn((indep)*alp); //Prob. function probit.
dpsest=pdfn((indep)*alp); //density. function probit.
psialp = (indep.*dpsest.*(inst-psest)./(psest.*(1-psest))); 
delpsi = dpsest.*indep;       
endif;
    

    Emu1alp=-meanc(((inst.*(dep-mu1))./(psest.^2)).*delpsi);
    Emu0alp=meanc((((1-inst).*(dep-mu0))./((1-psest).^2)).*delpsi);
    Em1alp=-meanc(((inst.*(treat-mu1d))./(psest.^2)).*delpsi);
    Em0alp=meanc((((1-inst).*(treat-mu0d))./((1-psest).^2)).*delpsi);
EH = -(1/N)*psialp'*psialp;

psimu1=(inst.*(dep-mu1)./psest); //v2
psimu0=((1-inst).*(dep-mu0)./(1-psest)); //v3
psim1=(inst.*(treat-mu1d)./psest); //v4
psim0=((1-inst).*(treat-mu0d)./(1-psest)); //v5

local nix1, nix2, nix3, Vnix;
nix1=-((1/gammah).*(Emu1alp-Emu0alp)-(deltah/gammah^2).*(Em1alp-Em0alp))'*inv(-EH)*((1/gammah).*(Emu1alp-Emu0alp)-(deltah/gammah^2).*(Em1alp-Em0alp));
nix2=meanc(((1/gammah).*(psimu1)-(deltah/gammah^2).*(psim1)).^2);
nix3=meanc(((1/gammah).*(psimu0)-(deltah/gammah^2).*(psim0)).^2);
Vnix=nix1+nix2+nix3;

local v_delta, v_gamma, cv_deltagamma, av_late, se_late;

  v_delta=-(Emu1alp-Emu0alp)'*inv(-EH)*(Emu1alp-Emu0alp)+meanc((inst.*(dep-mu1)./psest).^2)+meanc(((1-inst).*(dep-mu0)./(1-psest)).^2);
  v_gamma=-(Em1alp -Em0alp)'*inv(-EH)*(Em1alp-Em0alp)+meanc((inst.*(treat-mu1d)./psest).^2)+meanc(((1-inst).*(treat-mu0d)./(1-psest)).^2);
  cv_deltagamma=-(Emu1alp-Emu0alp)'*inv(-EH)*(Em1alp-Em0alp)+(1/N)*((inst.*(dep-mu1)./psest)'*(inst.*(treat-mu1d)./psest))+(1/N)*(((1-inst).*(dep-mu0)./(1-psest))'*((1-inst).*(treat-mu0d)./(1-psest)));
    av_late=(1/gammah)^2*v_delta+(tauh^2/gammah^2)*v_gamma-2*(tauh/gammah^2)*cv_deltagamma;
    se_late=sqrt((1/N)*av_late);
//sqrt((1/N)*Vnix),
retp(se_late);
endp;

