clear 
import com.comsol.model.util.*
ModelUtil.showProgress(true); 

% define inclusion parameters
% geometry length unit is mm
px=0.5; py=0.5; pz=2;  
rx=.10; ry=.10; rz=0.1;

% define model parameters
dT=10; 
I0=1e-2; 
len=3e-3; 
wd=1e-3; 


Nsph=10; % number of inclusions
Nbatch=2; % # of geometry realization 



rng('shuffle');
s=rng; 

PfBatch=zeros(Nbatch,1); 
SBatch=0*PfBatch;
sigAvBatch=0*PfBatch;
sigBatch=0*PfBatch; 
lamBatch=0*PfBatch; 
zTBatch=0*PfBatch; 
zTcmlBatch=0*PfBatch; 

ib=1; 
while ib<=Nbatch
    % load model
    model=mphload("thermoelectric_generator.noSpheres.vp04b.mph");


    model.param.set('L_element', [num2str(len),'[m]']);
    model.param.set('w_element', [num2str(wd),'[m]']);
    model.param.set('T0', [num2str(dT),'[degC]']);
    model.param.set('I0', [num2str(0),'[A]']);
    
    geom=model.component('comp1').geom('geom1'); 
    lblLst=cell(Nsph,1); 
    
    for isph=1:Nsph
        lbl=['e',num2str(isph)];

        % generate random locations of inclusions
        px=1.5*rx+(1-3*rx)*rand();
        py=1.5*ry+(1-3*ry)*rand();
        pz=0.1+1.5*rz+(2.8-3*rz)*rand();
    
        % put inclusions on a grid
        stp=1e2; 
        px=round(px*stp)/stp; 
        py=round(py*stp)/stp; 
        pz=round(pz*stp)/stp; 
    
        %insert inclusions into the geometry
        geom.create(lbl,'Ellipsoid'); 
        geom.feature(lbl).set('semiaxes',[rx,ry,rz]*1e-0); 
        geom.feature(lbl).set('pos', [px,py,pz]*1e-0); 
        geom.feature(lbl).set('rot', 30); 
        geom.feature(lbl).set('contributeto','csel1'); 
        lblLst{isph}=lbl; 
    end 
    geom.feature.create('uni1','Union'); 
    geom.feature('uni1').selection('input').set(lblLst);
    geom.feature('uni1').set('intbnd', 'off');
    geom.feature('uni1').set('contributeto','csel1'); 
    geom.run; 


    % run thermoelectric study 
    model.study('std1').run();
    S=abs(mphglobal(model,'ec.fp1.V0','dataset','dset1'))/dT; % extract Seebeck coefficient

    % run conductivity study
    model.param.set('I0', [num2str(I0),'[A]']);
    model.study('std2').run();

    V0=mphglobal(model,'ec.fp1.V0','dataset','dset2'); 
    sig=I0/V0*len/wd^2; % calculate electric conductivity
    q0=mphint2(model,'ht.ntflux','surface', 'selection', [3],'dataset','dset2'); 
    qCheck=mphint2(model,'ht.ntflux','surface', 'selection', [6],'dataset','dset2'); 

    % check for convergence
    if abs((q0-qCheck)/q0)<1e-2
        lam=q0*len/wd^2/dT; % thermal conductivity 
        Pf=S^2*sig; % power factor
        zf=Pf/lam; 
        zT=zf*273; %zT

        zTcmlBatch(ib)=mphmean(model,'tee1.Zmean*T','volume','dataset','dset1'); 
        sigAvBatch(ib)=mphmean(model,'ec.sigmazz','volume','dataset','dset2');
        SBatch(ib)=S; 
        PfBatch(ib)=Pf; 
        sigBatch(ib)=sig; 
        lamBatch(ib)=lam; 
        zTBatch(ib)=zT;
        disp(['advancing iteration, ib=',num2str(ib)]); 
        ib=ib+1;
    else 
        disp(['skipping iteration, ib=',num2str(ib)]); 
        mphsave(model,['./batchData/lamAnalysis/Nsph=',num2str(Nsph),'.ib=',num2str(ib),'.mph']);
        tmp=rand(); 
    end 

end 
save(['./Nsph=',num2str(Nsph),'.batchSize=',num2str(Nbatch),'.seed=',num2str(s.Seed),'.mat'], ...
    "SBatch", "PfBatch","sigBatch","lamBatch","zTBatch","zTcmlBatch","sigAvBatch"); 

