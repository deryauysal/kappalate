//asymptotic variance for normalized tau


proc(1)= asyvar_normkappalate(struct PV p, struct DS ds_struct); 
local indep, dep, treat, inst, N, pf ;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   inst = ds_struct[4].dataMatrix;  
   pf = ds_struct[5].Type; //probit or logit for the IPS 
   N = ds_struct[6].Type;

local alp, delta1h, delta0h, gamma1h, gamma0h, tauh;
     alp = pvUnpack(p, "alphahat");
     delta1h = pvUnpack(p, "delta1hat");
     gamma1h = pvUnpack(p, "gamma1hat");
     delta0h = pvUnpack(p, "delta0hat");
     gamma0h = pvUnpack(p, "gamma0hat");
     tauh = pvUnpack(p, "tauhat");
    
local delpsi, psest, dpsest, EH,Egamma1alp,Egamma0alp, Edelta1alp, Edelta0alp,psialp, psidelta1, psidelta0, psigamma1, psigamma0;

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



Edelta1alp=-meanc(((treat.*inst.*dep)./(psest.^2)).*delpsi)-meanc(((treat.*(1-inst).*dep)./((1-psest).^2)).*delpsi);
Edelta0alp=-meanc((((treat-1).*inst.*dep)./(psest.^2)).*delpsi)-meanc((((treat-1).*(1-inst).*dep)./((1-psest).^2)).*delpsi);
Egamma1alp=-meanc(((treat.*inst)./(psest.^2)).*delpsi)-meanc(((treat.*(1-inst))./((1-psest).^2)).*delpsi);
Egamma0alp=-meanc((((treat-1).*inst)./(psest.^2)).*delpsi)-meanc((((treat-1).*(1-inst))./((1-psest).^2)).*delpsi);

//psialp = (indep.*(inst-psest)); //Prob. function logit.
psidelta1 = (treat.*((inst-psest)./(psest.*(1-psest))).*dep-delta1h);
psigamma1 = (treat.*((inst-psest)./(psest.*(1-psest))) - gamma1h);
psidelta0 = ((1-treat).*(((1-inst)-(1-psest))./(psest.*(1-psest))).*dep- delta0h);
psigamma0 = ((1-treat).*(((1-inst)-(1-psest))./(psest.*(1-psest)))- gamma0h);
EH = -(1/N)*psialp'*psialp;
local nix11, nix1, nix2, nix3, Vnix;
nix11=Edelta1alp./gamma1h-Edelta0alp./gamma0h-(delta1h./gamma1h^2).*Egamma1alp+(delta0h./gamma0h^2).*Egamma0alp;
nix1=-nix11'*inv(-EH)*nix11;
nix2=meanc(((1/gamma1h).*(psidelta1)-(delta1h/gamma1h^2).*(psigamma1)).^2);
nix3=meanc(((1/gamma0h).*(psidelta0)-(delta0h/gamma0h^2).*(psigamma0)).^2);
Vnix=nix1+nix2+nix3;
/*
local v_delta, v_gamma, cv_deltagamma, av_late, se_late;

  v_delta=-(Emu1alp-Emu0alp)'*inv(-EH)*(Emu1alp-Emu0alp)+meanc((inst.*(dep-mu1)./psest).^2)+meanc(((1-inst).*(dep-mu0)./(1-psest)).^2);
  v_gamma=-(Em1alp -Em0alp)'*inv(-EH)*(Em1alp-Em0alp)+meanc((inst.*(treat-mu1d)./psest).^2)+meanc(((1-inst).*(treat-mu0d)./(1-psest)).^2);
  cv_deltagamma=-(Emu1alp-Emu0alp)'*inv(-EH)*(Em1alp-Em0alp)+(1/N)*((inst.*(dep-mu1)./psest)'*(inst.*(treat-mu1d)./psest))+(1/N)*(((1-inst).*(dep-mu0)./(1-psest))'*((1-inst).*(treat-mu0d)./(1-psest)));
    av_late=(1/gammah)^2*v_delta+(tauh^2/gammah^2)*v_gamma-2*(tauh/gammah^2)*cv_deltagamma;
    se_late=sqrt((1/N)*av_late);*/
//se_late,
retp(sqrt((1/N)*Vnix));
endp;

