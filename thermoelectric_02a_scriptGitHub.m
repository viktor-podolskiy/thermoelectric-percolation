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


dx=0.025; % spheres will be placed on the grid with this mesh step, in mm
Nx=floor((wd*1e3-3*rx)/dx); 
Nz=floor((len*1e3-3*rz-0.2)/dx); 
xArr=1.5*rx+dx*(1:Nx); 
zArr=0.1+1.5*rz+dx*(1:Nz); 
[xS,yS,zS]=meshgrid(xArr,xArr,zArr); 

df=0.01; % finer mesh will be used for analyzing the structure of conductive clusters
[xF,yF,zF]=meshgrid((0:df:wd*1e3),(0:df:wd*1e3),(0.1:df:len*1e3-0.1)); 


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
numClust=0*PfBatch; 
zlenClust=0*PfBatch; 
zlenTot=0*PfBatch; 

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

        % mesh-based implementation
        gridIndex = randi(numel(xS));
        px = xS(gridIndex);
        py = yS(gridIndex);
        pz = zS(gridIndex);

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

        % analyze cluster structure
        xyzArr=[xF(:).'; yF(:).'; zF(:).']; 
        sig3DData=mphinterp(model,'ec.sigmazz','coord',xyzArr,'dataset','dset2'); 
        sig3DData = reshape(sig3DData, size(xF));
        thr=max(sig3DData(:))*0.8; 
        sig3DData(sig3DData<thr)=0; 
        sig3DData(sig3DData>=thr)=1; 

        tmp = bwconncomp(sig3DData); 
        numClust(ib)=tmp.NumObjects; 

        bbData=regionprops(tmp,"BoundingBox");
        bbMax=0; 
        for ibb=1:length(bbData)
            bbox = bbData(ibb).BoundingBox;
            if bbMax<bbox(6)
                bbMax=bbox(6); 
            end
        end
        zlenClust(ib)=bbMax; 
        zlenTot(ib)=size(xF,3); 

        save(['./Nsph=',num2str(Nsph),'.batchSize=',num2str(Nbatch),'.seed=',num2str(s.Seed),'.mat'], ...
            "SBatch", "PfBatch","sigBatch","lamBatch","zTBatch","zTcmlBatch","sigAvBatch","numClust","zlenTot","zlenClust"); 
        disp(['advancing iteration, ib=',num2str(ib), '; numClust=',num2str(tmp.NumObjects)]); 
        ib=ib+1;
    else 
        disp(['skipping iteration, ib=',num2str(ib)]); 
        mphsave(model,['./batchData/lamAnalysis/Nsph=',num2str(Nsph),'.ib=',num2str(ib),'.mph']);
        tmp=rand(); 
    end 

end 

