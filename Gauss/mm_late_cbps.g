

struct OutputsLATE {
   matrix late;
   matrix latese;
   matrix pars; // ips coeff, to means in numerator tow means in denominator
   matrix parsse;
   matrix retcode;
};



proc(1) =  mm_late(&f1, &f2,  struct DS ds_struct); //&f3,
    library dc;
local indep, dep, treat, inst, pf, N, v_mest, av_late, se_late;
   indep = ds_struct[1].dataMatrix; 
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   inst = ds_struct[4].dataMatrix; 
   pf = ds_struct[5].Type;
   N = ds_struct[6].Type;

 //Data, mean spec. linear: mf==1, logit: 2 poisson: 3
    //local alp_start, estcoef,gamhat, deltahat, latehat;
    local f1:proc; //one of the 5 
    local f2:proc; // procedure for crossproduct of the moments
    //local f3: proc; //procedure for asymm. variance
    
    //starting values
    local alp_start,psest, mu1s, mu0s, m1s, m0s, deltas, gammas,taus, delta1s, delta0s, gamma1s, gamma0s;
    
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
    if pf==1;
        dcout1 = binaryLogit(dcCt);
    elseif pf==2;
        dcout1 = binaryProbit(dcCt);
    endif;
    alp_start=pvGetParVector(dcout1.par);
if pf==1;//if logit
psest=exp((indep)*alp_start)./(1+exp((indep)*alp_start)); //Prob. function logit.
elseif pf==2;//if probit
psest=cdfn((indep)*alp_start); //Prob. function probit.
endif;    

    mu1s=meanc(inst.*(dep)./psest)/meanc((inst)./psest);     //   meanc((inst./psest).*treat); //meanc(dep.*inst);
    mu0s=meanc((1-inst).*(dep)./(1-psest))/meanc((1-inst)./(1-psest));   //meanc(((1-inst)./(1-psest)).*dep);//meanc(dep.*(1-inst))
    m1s =meanc(inst.*(treat)./psest)/meanc((inst)./psest); //  meanc((inst./psest).*treat); //meanc(treat.*inst);
    m0s =meanc((1-inst).*(treat)./(1-psest))/meanc((1-inst)./(1-psest));  //meanc(((1-inst)./(1-psest)).*treat); //meanc(treat.*(1-inst));

    if &f1==&normlate_mom_cbps;
    deltas=mu1s-mu0s;
    else;
    deltas=meanc((inst./psest).*dep-((1-inst)./(1-psest)).*dep);
    endif;

    if &f1==&kappa1late_mom_cbps;
    gammas=meanc((inst./psest).*treat-((1-inst)./(1-psest)).*treat);
    elseif  &f1==&kappa0late_mom_cbps;  
    gammas=meanc((inst./psest).*(treat-1)-((1-inst)./(1-psest)).*(treat-1));  
    elseif  &f1==&kappalate_mom_cbps;  
    gammas=meanc(1-(treat.*(1-inst))./(1-psest)-((1-treat).*inst)./psest); 
    elseif &f1==&normlate_mom_cbps;   
    gammas=  m1s-m0s;
    endif;    
    delta1s=meanc(treat.*((inst-psest)./(psest.*(1-psest))).*dep);
    delta0s=meanc((1-treat).*(((1-inst)-(1-psest))./(psest.*(1-psest))).*dep);
    gamma1s=meanc(treat.*((inst-psest)./(psest.*(1-psest))));
    gamma0s=meanc((1-treat).*(((1-inst)-(1-psest))./(psest.*(1-psest))));
    if &f1==&normkappalate_mom_cbps;
    taus= (delta1s/gamma1s-delta0s/gamma0s);   
    else;  
    taus= deltas/gammas;   
    endif;    
    
struct PV par;
par=pvCreate();

   //Regression
    if &f1==&kappa1late_mom_cbps;
par= pvPack(par,alp_start, "alphahat");
par= pvPack(par, deltas, "deltahat"); //num
par= pvPack(par, gammas, "gammahat"); //denom
par= pvPack(par, taus, "tauhat"); //denom      
        
    elseif &f1== &kappa0late_mom_cbps;
par= pvPack(par,alp_start, "alphahat");
par= pvPack(par, deltas, "deltahat"); //num
par= pvPack(par, gammas, "gammahat"); //denom
par= pvPack(par, taus, "tauhat"); //denom              
        
    elseif &f1== &kappalate_mom_cbps;
par= pvPack(par,alp_start, "alphahat");
par= pvPack(par, deltas, "deltahat"); //num
par= pvPack(par, gammas, "gammahat"); //denom
par= pvPack(par, taus, "tauhat"); //denom      
        
    elseif &f1== &normkappalate_mom_cbps;    
par= pvPack(par,alp_start, "alphahat");
par= pvPack(par, delta1s, "delta1hat"); //num
par= pvPack(par, gamma1s, "gamma1hat"); //denom
par= pvPack(par, delta0s, "delta0hat"); //num
par= pvPack(par, gamma0s, "gamma0hat"); //denom
par= pvPack(par, taus, "tauhat"); //denom      

   
       
    elseif &f1== &normlate_mom_cbps;
par= pvPack(par,alp_start, "alphahat");
par= pvPack(par, mu1s, "mu1hat"); 
par= pvPack(par, mu0s, "mu0hat");   
par= pvPack(par, m1s, "m1hat"); 
par= pvPack(par, m0s, "m0hat");        
par= pvPack(par, deltas, "deltahat"); //num
par= pvPack(par, gammas, "gammahat"); //denom
par= pvPack(par, taus, "tauhat"); //denom           
    endif;    
local latehat, estcoef;

// Declare control structure and fill with defaults
struct eqSolvemtControl c;
c = eqSolvemtControlCreate();
c.maxIters=10000;    
c.output = 0;
// Declare output structure to hold results
struct eqSolvemtOut out;
// Solve the system of equations
out =  eqSolvemt(&f1, par, ds_struct, c);

    latehat=pvUnpack(out.par, "tauhat");  //late estimate
    estcoef=pvGetParVector(out.par); //all coeffecient estimates

    
    v_mest=(1/N)*inv(gradMT(&f1,out.par, ds_struct))*f2(out.par, ds_struct)*inv(gradMT(&f1,out.par, ds_struct))';
    /*v_delta=v_mest[rows(estcoef)-1,rows(estcoef)-1];
    v_gamma=v_mest[rows(estcoef),rows(estcoef)];
    cv_deltagamma=v_mest[rows(estcoef),rows(estcoef)-1];
    av_late=(1/gamhat)^2*v_delta+(latehat^2/gamhat^2)*v_gamma-2*(latehat/gamhat^2)*cv_deltagamma;
    */
    av_late= v_mest[rows(estcoef),rows(estcoef)];
    se_late=sqrt(av_late);
    //Define structure

   struct OutputsLATE outlate;
   outlate.late  = latehat;
   outlate.latese  = se_late;
   outlate.pars=  estcoef;
   outlate.parsse = sqrt(diag(v_mest));
   outlate.retcode = out.retcode;
   retp(outlate);
 
 
    //retp(ate, se_mest, estcoef, sqrt(diag(v_mest)), out.par,out.retcode);


endp;

proc kappa1late_crossmom_cbps(struct PV p, struct DS ds_struct);
local psi1, alp, gammah, deltah, tauh, psi2, psi3, psi4;
     alp = pvUnpack(p, "alphahat");
     gammah = pvUnpack(p, "gammahat");
     deltah = pvUnpack(p, "deltahat");
     tauh = pvUnpack(p, "tauhat");
    
    local indep, dep, treat, inst, pf, N;
   indep = ds_struct[1].dataMatrix; 
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   inst = ds_struct[4].dataMatrix; 
   pf = ds_struct[5].Type;
   N = ds_struct[6].Type;
local   psest, dpsest, v11,v12,v13,v14, v21,v22, v23, v24, v31, v32, v33, v34, v41, v42, v43, v44;

if pf==1;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = (indep.*(inst-psest)./(psest.*(1-psest))); //Prob. function logit.

endif;

    
psi2=((inst./psest).*dep-((1-inst)./(1-psest)).*dep)-deltah;
psi3=((inst./psest).*treat-((1-inst)./(1-psest)).*treat)-gammah;
psi4=ones(n,1).*(deltah/gammah-tauh);


v11=(1/N)*psi1'*psi1;
v12=(1/N)*psi1'*psi2;
v13=(1/N)*psi1'*psi3;
v14=(1/N)*psi1'*psi4;

v21=(1/N)*psi2'*psi1;
v22=(1/N)*psi2'*psi2;
v23=(1/N)*psi2'*psi3;
v24=(1/N)*psi2'*psi4;

v31=(1/N)*psi3'*psi1;
v32=(1/N)*psi3'*psi2;
v33=(1/N)*psi3'*psi3;
v34=(1/N)*psi3'*psi4;

v41=(1/N)*psi4'*psi1;
v42=(1/N)*psi4'*psi2;
v43=(1/N)*psi4'*psi3;
v44=(1/N)*psi4'*psi4;

retp ((v11~v12~v13~v14)|(v21~v22~v23~v24)|(v31~v32~v33~v34)|(v41~v42~v43~v44));

endp;


proc kappa0late_crossmom_cbps(struct PV p, struct DS ds_struct);
local psi1, alp, gammah, deltah, tauh, psi2, psi3, psi4;
     alp = pvUnpack(p, "alphahat");
     gammah = pvUnpack(p, "gammahat");
     deltah = pvUnpack(p, "deltahat");
     tauh = pvUnpack(p, "tauhat");

    local indep, dep, treat, inst, pf, N;
   indep = ds_struct[1].dataMatrix; 
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   inst = ds_struct[4].dataMatrix; 
   pf = ds_struct[5].Type;
   N = ds_struct[6].Type;

local   psest, dpsest, v11,v12,v13,v14, v21,v22, v23, v24, v31, v32, v33, v34, v41, v42, v43, v44;

if pf==1;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = (indep.*(inst-psest)./(psest.*(1-psest))); //Prob. function logit.
endif;


psi2=((inst./psest).*dep-((1-inst)./(1-psest)).*dep)-deltah;
//psi3=(((1-inst)./(1-psest))-(inst./psest)+(inst./psest).*treat-((1-inst)./(1-psest)).*treat)-gammah;
psi3=((inst./psest).*(treat-1)-((1-inst)./(1-psest)).*(treat-1))-gammah;
psi4=ones(n,1).*(deltah/gammah-tauh);


v11=(1/N)*psi1'*psi1;
v12=(1/N)*psi1'*psi2;
v13=(1/N)*psi1'*psi3;
v14=(1/N)*psi1'*psi4;

v21=(1/N)*psi2'*psi1;
v22=(1/N)*psi2'*psi2;
v23=(1/N)*psi2'*psi3;
v24=(1/N)*psi2'*psi4;

v31=(1/N)*psi3'*psi1;
v32=(1/N)*psi3'*psi2;
v33=(1/N)*psi3'*psi3;
v34=(1/N)*psi3'*psi4;

v41=(1/N)*psi4'*psi1;
v42=(1/N)*psi4'*psi2;
v43=(1/N)*psi4'*psi3;
v44=(1/N)*psi4'*psi4;

retp ((v11~v12~v13~v14)|(v21~v22~v23~v24)|(v31~v32~v33~v34)|(v41~v42~v43~v44));

endp;

       


proc kappalate_crossmom_cbps(struct PV p, struct DS ds_struct);
local psi1, alp, gammah, deltah, tauh, psi2, psi3, psi4;
     alp = pvUnpack(p, "alphahat");
     gammah = pvUnpack(p, "gammahat");
     deltah = pvUnpack(p, "deltahat");
     tauh = pvUnpack(p, "tauhat");

    local indep, dep, treat, inst, pf, N;
   indep = ds_struct[1].dataMatrix; 
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   inst = ds_struct[4].dataMatrix; 
   pf = ds_struct[5].Type;
   N = ds_struct[6].Type;

local   dpsest, psest, v11,v12,v13,v14, v21,v22, v23, v24, v31, v32, v33, v34, v41, v42, v43, v44;

if pf==1;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = (indep.*(inst-psest)./(psest.*(1-psest))); //Prob. function logit.
endif;


psi2=((inst./psest).*dep-((1-inst)./(1-psest)).*dep)-deltah;
//psi3=(1-(inst./psest)+(inst./psest).*treat-((1-inst)./(1-psest)).*treat)-gammah;
psi3=1-(treat.*(1-inst))./(1-psest)-((1-treat).*inst)./psest-gammah;
psi4=ones(n,1).*(deltah/gammah-tauh);


v11=(1/N)*psi1'*psi1;
v12=(1/N)*psi1'*psi2;
v13=(1/N)*psi1'*psi3;
v14=(1/N)*psi1'*psi4;

v21=(1/N)*psi2'*psi1;
v22=(1/N)*psi2'*psi2;
v23=(1/N)*psi2'*psi3;
v24=(1/N)*psi2'*psi4;

v31=(1/N)*psi3'*psi1;
v32=(1/N)*psi3'*psi2;
v33=(1/N)*psi3'*psi3;
v34=(1/N)*psi3'*psi4;

v41=(1/N)*psi4'*psi1;
v42=(1/N)*psi4'*psi2;
v43=(1/N)*psi4'*psi3;
v44=(1/N)*psi4'*psi4;

retp ((v11~v12~v13~v14)|(v21~v22~v23~v24)|(v31~v32~v33~v34)|(v41~v42~v43~v44));


endp;



   
proc normkappalate_crossmom_cbps(struct PV p, struct DS ds_struct);

   local alp, psi1,psi2,psi3, psi4,psi5,psi6, delta1h, delta0h, gamma1h, gamma0h, tauh;
     alp = pvUnpack(p, "alphahat");
     delta1h = pvUnpack(p, "delta1hat");
     gamma1h = pvUnpack(p, "gamma1hat");
     delta0h = pvUnpack(p, "delta0hat");
     gamma0h = pvUnpack(p, "gamma0hat");
     tauh = pvUnpack(p, "tauhat");
    
    local indep, dep, treat, inst, pf, N;
   indep = ds_struct[1].dataMatrix; 
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   inst = ds_struct[4].dataMatrix; 
   pf = ds_struct[5].Type;
   N = ds_struct[6].Type;

local   psest,dpsest, v11,v12,v13,v14, v15, v16, v21,v22, v23, v24, v25, v26, v31, v32, v33, v34, v35, v36;
local   v41,v42,v43,v44, v45, v46, v51,v52, v53, v54, v55, v56, v61, v62, v63, v64, v65, v66;


if pf==1;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = (indep.*(inst-psest)./(psest.*(1-psest))); //Prob. function logit.
endif;


psi2 = (treat.*((inst-psest)./(psest.*(1-psest))).*dep-delta1h);
psi3 = (treat.*((inst-psest)./(psest.*(1-psest))) - gamma1h);
psi4 = ((1-treat).*(((1-inst)-(1-psest))./(psest.*(1-psest))).*dep- delta0h);
psi5 = ((1-treat).*(((1-inst)-(1-psest))./(psest.*(1-psest)))- gamma0h);
psi6 = ones(n,1).*(delta1h/gamma1h-delta0h/gamma0h-tauh);    


v11=(1/N)*psi1'*psi1;
v12=(1/N)*psi1'*psi2;
v13=(1/N)*psi1'*psi3;
v14=(1/N)*psi1'*psi4;
v15=(1/N)*psi1'*psi5;
v16=(1/N)*psi1'*psi6;


v21=(1/N)*psi2'*psi1;
v22=(1/N)*psi2'*psi2;
v23=(1/N)*psi2'*psi3;
v24=(1/N)*psi2'*psi4;
v25=(1/N)*psi2'*psi5;
v26=(1/N)*psi2'*psi6;

v31=(1/N)*psi3'*psi1;
v32=(1/N)*psi3'*psi2;
v33=(1/N)*psi3'*psi3;
v34=(1/N)*psi3'*psi4;
v35=(1/N)*psi3'*psi5;
v36=(1/N)*psi3'*psi6;

v41=(1/N)*psi4'*psi1;
v42=(1/N)*psi4'*psi2;
v43=(1/N)*psi4'*psi3;
v44=(1/N)*psi4'*psi4;
v45=(1/N)*psi4'*psi5;
v46=(1/N)*psi4'*psi6;

v51=(1/N)*psi5'*psi1;
v52=(1/N)*psi5'*psi2;
v53=(1/N)*psi5'*psi3;
v54=(1/N)*psi5'*psi4;
v55=(1/N)*psi5'*psi5;
v56=(1/N)*psi5'*psi6;

v61=(1/N)*psi6'*psi1;
v62=(1/N)*psi6'*psi2;
v63=(1/N)*psi6'*psi3;
v64=(1/N)*psi6'*psi4;
v65=(1/N)*psi6'*psi5;
v66=(1/N)*psi6'*psi6;

retp ((v11~v12~v13~v14~v15~v16)|(v21~v22~v23~v24~v25~v26)|(v31~v32~v33~v34~v35~v36)|(v41~v42~v43~v44~v45~v46)|(v51~v52~v53~v54~v55~v56)|(v61~v62~v63~v64~v65~v66));

endp;



proc normlate_crossmom_cbps(struct PV p, struct DS ds_struct);
local alp, psi1,psi2,psi3, psi4,psi5,psi6,psi7,psi8;
local psest, dpsest , mu1, mu0, mu1d, mu0d, gammah, deltah, tauh;
     alp = pvUnpack(p, "alphahat");
     mu1 = pvUnpack(p, "mu1hat");
     mu0 = pvUnpack(p, "mu0hat");
     mu1d = pvUnpack(p, "m1hat");
     mu0d = pvUnpack(p, "m0hat");
     deltah = pvUnpack(p, "deltahat");
     gammah = pvUnpack(p, "gammahat");
     tauh = pvUnpack(p, "tauhat");
    
    local indep, dep, treat, inst, N, pf;
   indep = ds_struct[1].dataMatrix; 
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   inst = ds_struct[4].dataMatrix; 
   pf = ds_struct[5].Type;
   N = ds_struct[6].Type;
    
local   v11,v12,v13,v14, v15, v16, v17, v18, v21,v22, v23, v24, v25, v26, v27, v28, v31, v32, v33, v34, v35, v36, v37, v38;
local   v41,v42,v43,v44, v45, v46, v47, v48, v51,v52, v53, v54, v55, v56, v57, v58, v61, v62, v63, v64, v65, v66, v67, v68;
local  v71, v72, v73, v74, v75, v76, v77, v78;
local  v81, v82, v83, v84, v85, v86, v87, v88; 
local v1, v2, v3, v4, v5 ,v6, v7, v8;

if pf==1;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = (indep.*(inst-psest)./(psest.*(1-psest))); //Prob. function logit.
endif;


psi2=(inst.*(dep-mu1)./psest);
psi3=((1-inst).*(dep-mu0)./(1-psest));
psi4=(inst.*(treat-mu1d)./psest);
psi5=((1-inst).*(treat-mu0d)./(1-psest));
psi6=ones(n,1).*(mu1-mu0-deltah);
psi7=ones(n,1).*(mu1d-mu0d-gammah);
psi8 = ones(n,1).*(deltah/gammah-tauh);   


v11=(1/N)*psi1'*psi1;
v12=(1/N)*psi1'*psi2;
v13=(1/N)*psi1'*psi3;
v14=(1/N)*psi1'*psi4;
v15=(1/N)*psi1'*psi5;
v16=(1/N)*psi1'*psi6;
v17=(1/N)*psi1'*psi7;
v18=(1/N)*psi1'*psi8;

v21=(1/N)*psi2'*psi1;
v22=(1/N)*psi2'*psi2;
v23=(1/N)*psi2'*psi3;
v24=(1/N)*psi2'*psi4;
v25=(1/N)*psi2'*psi5;
v26=(1/N)*psi2'*psi6;
v27=(1/N)*psi2'*psi7;
v28=(1/N)*psi2'*psi8;

v31=(1/N)*psi3'*psi1;
v32=(1/N)*psi3'*psi2;
v33=(1/N)*psi3'*psi3;
v34=(1/N)*psi3'*psi4;
v35=(1/N)*psi3'*psi5;
v36=(1/N)*psi3'*psi6;
v37=(1/N)*psi3'*psi7;
v38=(1/N)*psi3'*psi8;

v41=(1/N)*psi4'*psi1;
v42=(1/N)*psi4'*psi2;
v43=(1/N)*psi4'*psi3;
v44=(1/N)*psi4'*psi4;
v45=(1/N)*psi4'*psi5;
v46=(1/N)*psi4'*psi6;
v47=(1/N)*psi4'*psi7;
v48=(1/N)*psi4'*psi8;

v51=(1/N)*psi5'*psi1;
v52=(1/N)*psi5'*psi2;
v53=(1/N)*psi5'*psi3;
v54=(1/N)*psi5'*psi4;
v55=(1/N)*psi5'*psi5;
v56=(1/N)*psi5'*psi6;
v57=(1/N)*psi5'*psi7;
v58=(1/N)*psi5'*psi8;

v61=(1/N)*psi6'*psi1;
v62=(1/N)*psi6'*psi2;
v63=(1/N)*psi6'*psi3;
v64=(1/N)*psi6'*psi4;
v65=(1/N)*psi6'*psi5;
v66=(1/N)*psi6'*psi6;
v67=(1/N)*psi6'*psi7;
v68=(1/N)*psi6'*psi8;

v71=(1/N)*psi7'*psi1;
v72=(1/N)*psi7'*psi2;
v73=(1/N)*psi7'*psi3;
v74=(1/N)*psi7'*psi4;
v75=(1/N)*psi7'*psi5;
v76=(1/N)*psi7'*psi6;
v77=(1/N)*psi7'*psi7;
v78=(1/N)*psi7'*psi8;

v81=(1/N)*psi8'*psi1;
v82=(1/N)*psi8'*psi2;
v83=(1/N)*psi8'*psi3;
v84=(1/N)*psi8'*psi4;
v85=(1/N)*psi8'*psi5;
v86=(1/N)*psi8'*psi6;
v87=(1/N)*psi8'*psi7;
v88=(1/N)*psi8'*psi8;

v1= (v11~v12~v13~v14~v15~v16~v17~v18);
v2= (v21~v22~v23~v24~v25~v26~v27~v28);
v3= (v31~v32~v33~v34~v35~v36~v37~v38);
v4= (v41~v42~v43~v44~v45~v46~v47~v48);
v5= (v51~v52~v53~v54~v55~v56~v57~v58);
v6= (v61~v62~v63~v64~v65~v66~v67~v68);
v7= (v71~v72~v73~v74~v75~v76~v77~v78);
v8= (v81~v82~v83~v84~v85~v86~v87~v88);

retp (v1|v2|v3|v4|v5|v6|v7|v8);


endp;
