%function [] = ClathrinMediatedEndocytosis(dirRoot,dirMod,varargin)
%--------------------------------------------------------------------------
        % ClathrinMediatedEndocytosis performs the dynamics of CALM,
        % Clathrin and membrane during endocytosis
        % input: 
        % dirRoot - root directory for data storage
        % dirMod - directory of @ModMembrane and other required objects
        % optional:
        % see variable arguments
        %   See also tetherPull, membraneFussion, exampleLamellipodia,
        %   exampleEndocytosis, redBloodCell
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------        
% ip = inputParser;
% ip.CaseSensitive = false;
% ip.addRequired('dirRoot', @(x) ischar(x));
% ip.addRequired('dirMod', @(x) ischar(x));
% ip.addParameter('nStep', 20000, @isnumeric); %simulation steps for membrane relaxation
% ip.addParameter('nRep', 12, @isnumeric); %repeat number
% ip.addParameter('xyzLim', [-20 20;-20 20;-20 20], @isnumeric); %mesh range
% ip.addParameter('xyzLimPlot', [-10 10;-10 10;-15 15], @isnumeric); %mesh range
% ip.addParameter('viewAng', [0 0], @isnumeric); %figure view angle
% ip.addParameter('idInit', [], @isnumeric); %which initial condition to use
% ip.addParameter('Tcutoff', 0.00000005, @isnumeric); %time scale of plotting
% ip.parse(dirRoot,dirMod,varargin{:});
% dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/CME2';
% dirMod='/home2/s171152/codes/matlab/mine/git/DanuserLab/biophysicsmodels';addpath(genpath(dirMod));
% nStep=2000;nRep=12;xyzLim=[-20 20;-20 20;-20 20];viewAng=[0 0];xyzLimPlot=[-10 10;-10 10;-15 15];idInit=1;
%--------------------------------------------------------------------------
% nStep=ip.Results.nStep;
% nRep=ip.Results.nRep;
% xyzLim=ip.Results.xyzLim;
% viewAng=ip.Results.viewAng;
% xyzLimPlot=ip.Results.xyzLimPlot;
% idInit=ip.Results.idInit;
% Tcutoff=ip.Results.Tcutoff;%0.00000005
%--------------------------------------------------------------------------
%% initial parameters
dirRoot = '/endosome/work/bioinformatics/s171152/data/CMEmodel';
dirRootInitMemCALM=[dirRoot filesep 'MemCALM'];
dirRootSt1=[dirRoot filesep 'St1'];
dirRootSt2=[dirRoot filesep 'St2'];
dirRootSt3=[dirRoot filesep 'St3'];
mkdir(dirRoot); mkdir(dirRootInitMemCALM); mkdir(dirRootSt1); mkdir(dirRootSt2); mkdir(dirRootSt3);
dirMod='/home2/s171152/codes/matlab/mine/git/DanuserLab/biophysicsmodels';addpath(genpath(dirMod));
xyzLim=[-20 20;-20 20;-20 20];viewAng=[0 0];xyzLimPlot=[-10 10;-10 10;-15 15];idInit=1;Tcutoff=0.00000005;
%%
% Initialization
%--------------------------------------------------------------------------
%setup parameters and variables
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
m=ModMembrane(true,4,0,'unit',u);
s=ModSubstrate(0.01,13);
c=ModClathrin('unit',u);
p=ModFreeParticle(0,'unit',u);
p.prop = {'Particle','needFollow'};

%assemble the model with the membrane object
int_info=struct('TypName',{'TypForce';'TypChemistry'},'IntList',{{{c,p},{m,s}};{{}}});
M = model(int_info,dirMod,u,xyzLim,m,s,c,p);
% computation of inital membrane with CALM
%--------------------------------------------------------------------------
% k_V=[1 8];
% k_A=[1 2 8];
P=[1];
kA=[8];
rec=ComRecord([dirRootInitMemCALM filesep 'data'],[dirRootInitMemCALM filesep 'fig'],[dirRootInitMemCALM filesep 'init'],...
              'nRep',4,'deleteFiles',true,...
              'paraAlt',struct('P',P,'kA',kA,'iCALM',[0 1 2]));
nPar=numel(rec.dir_all);
%--------------------------------------------------------------------------
[minZ,idTem]=min(M.mod{M.i_mod.ModMembrane}.var.coord(:,3));
% set CALM position
r_ctr=M.mod{M.i_mod.ModMembrane}.var.coord(idTem,:);
[~,id_tem]=sort(sum((m.var.coord-r_ctr).^2,2));
ScoordLow=M.mod{M.i_mod.ModMembrane}.var.coord(id_tem(1:19),:);%20
ScoordLow(:,:)=ScoordLow(:,:)*2;
ScoordLow(:,3)=-4*ScoordLow(:,3);
ScoordLow(:,3)=ScoordLow(:,3)-max(ScoordLow(:,3));
ScoordLow(:,3)=minZ+ScoordLow(:,3)+3;%+3.5
% ScoordLow=[ScoordLow; [ScoordLow(:,1),ScoordLow(:,2),-ScoordLow(:,3)]];

ScoordHigh=ScoordLow;
ScoordHigh(:,3)=2*ScoordHigh(:,3);
ScoordHigh(:,3)=ScoordHigh(:,3)-max(ScoordHigh(:,3));
ScoordHigh(:,3)=minZ+ScoordHigh(:,3)+5;

ScoordFlat=ScoordLow;
ScoordFlat(:,3)=min(ScoordLow(:,3));


M.mod{M.i_mod.ModSubstrate}.var.n_coord=size(ScoordLow,1);
% adjust membrane parameters
%--------------------------------------------------------------------------
M.mod{M.i_mod.ModMembrane}.pm.A0=4*pi*abs(minZ)^2*0.9;
M.mod{M.i_mod.ModMembrane}.pm.V0=4/3*pi*abs(minZ)^3*1.1;
M.mod{M.i_mod.ModMembrane}.pm.k_V=0;
M.mod{M.i_mod.ModMembrane}.pm.mu=500000; %500000
M.mod{M.i_mod.ModMembrane}.pm.Vdh.V0=0.1; %0.2
M.mod{M.i_mod.ModMembrane}.pm.f_const_std_std=0.0005;
M.mod{M.i_mod.ModMembrane}.pm.DlocRelax=20;
%--------------------------------------------------------------------------
%setup parallel variables to support parfor loop    
Mpar=cell(nPar,1);
dirPar=cell(nPar,1);
for iPar=1:nPar
    Mpar{iPar}=M;
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.k_V=rec.dir_all{iPar}.paraVal(1);
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.k_A=rec.dir_all{iPar}.paraVal(2);
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.P=rec.dir_all{iPar}.paraVal(1);
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.k_A=rec.dir_all{iPar}.paraVal(2);
    Mpar{iPar}.TypForce.pm.k_ModMembrane_ModSubstrate=10;  
    if rec.dir_all{iPar}.paraVal(3)==0
        Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModSubstrate}.var.coord=ScoordFlat;
    elseif rec.dir_all{iPar}.paraVal(3)==1
        Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModSubstrate}.var.coord=ScoordLow;
    else
        Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModSubstrate}.var.coord=ScoordHigh;
    end
    dirPar{iPar}=rec.dir_all{iPar};
end
%--------------------------------------------------------------------------
%%
%Computation
parfor iPar=1:nPar
    delete([dirPar{iPar}.dir_data filesep '*.*']);
    delete([dirPar{iPar}.dir_fig filesep '*.*']);
    i_loop_done=0;
    fig=ComPlot.PlotMod(Mpar{iPar}, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
    ComRecord.addFrame(dirPar{iPar},Mpar{iPar},fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false); 
end

nStep=10000;
for iPar=1:nPar
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModFreeParticle}.prop={'Particle'};
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.mu=500000;
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.DlocRelax=20;
end
%%---------------------------------------------------------------------
fprintf('computing CME morphology with CALM...\n');
parfor iPar=1:nPar
% for iPar=1:1
    M=Mpar{iPar};
    %setup dynamics computation, indicating simulation steps, type of
    %forces included, whether adaptive time step is needed and so forth
    dyn=dynamics(nStep);
    Fname={'ModMembrane';'ModMembrane_ModSubstrate'};
%     Fname={'ModClathrin';'ModMembrane';'ModMembrane_ModSubstrate';'ModClathrin_ModFreeParticle'};
    [dyn] = Preparation(dyn,M,Fname);
    i_loop_done=0;
    while (i_loop_done<nStep)
        [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',Tcutoff,'plot_or_not',false,'sMax', 0.1);
        i_loop_done=i_loop_done+CutOff.nS;
        %plot and save results
        fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
        ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false); 
        fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
    end 
Mpar{iPar}=M;
%makeMovie(rec,'DelayTime', 0.25);
fprintf('finishded at i_par %d \n',iPar);
end
%%
dirRootSt1Tem=[dirRootSt1 filesep num2str(1)];
rec=ComRecord([dirRootSt1Tem filesep 'data'],[dirRootSt1Tem filesep 'fig'],[dirRootSt1Tem filesep 'init'],...
              'nRep',12,'deleteFiles',true);
nPar=numel(rec.dir_all);
Mpar=cell(nPar,1);
dirPar=cell(nPar,1);
S=load([dirRootSt1Tem filesep 'FlatInit.mat']);
for iPar=1:4         
    Mpar{iPar}=S.M;
    dirPar{iPar}=rec.dir_all{iPar};
end
S=load([dirRootSt1Tem filesep 'LowInit.mat']);
for iPar=5:8
    Mpar{iPar}=S.M;
    dirPar{iPar}=rec.dir_all{iPar};
end
S=load([dirRootSt1Tem filesep 'HighInit.mat']);
for iPar=9:12
    Mpar{iPar}=S.M;
    dirPar{iPar}=rec.dir_all{iPar};
end
%--------------------------------------------------------------------------
CtrPar=cell(nPar,1);
nCtr=1;
for iPar=1:nPar
    R=Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.var.coord;
    idDown=R(:,3)<0;
    R=R(idDown,:);
    [~,idSort]=sort(vecnorm(R(:,1:2),2,2));
    Zctr=abs(R(idSort(1),3));
    
    CtrPar{iPar}=zeros(nCtr,3);
    CtrPar{iPar}(1,:)=[0,0,-(Zctr-3)];
%     CtrPar{iPar}(2,:)=[0,0,(Zctr-3)];
end
fprintf('initializing 1st Clathrin...\n');
parfor iPar=1:nPar
% for iPar=1:1 
    M=Mpar{iPar};
    for iCtr=1:nCtr
            if iCtr==1
                update=false;
            else
                update=true;
            end
    changed=false;
    while (changed==false)
        [M,changed] = init(M.mod{M.i_mod.ModClathrin},M,CtrPar{iPar}(iCtr,:),'update',update);
    end
    end
    Mpar{iPar}=M;
end
name={'ModClathrin_ModFreeParticle'};
parfor iPar=1:nPar
% for iPar=1:1  
    M=Mpar{iPar};
    [M] = locDyn(M.mod{M.i_mod.ModClathrin},M,name,1,(1:3),'sMax',0.01,'update',false,'initOtherInfo',true,'addOtherInfo',true);
    Mpar{iPar}=M;
    i_loop_done=0;
    fig=ComPlot.PlotMod(Mpar{iPar}, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
    ComRecord.addFrame(dirPar{iPar},Mpar{iPar},fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false); 
end

% save([dirRoot filesep 'stg1(1stClathrin).mat']);
%%---------------------------------------------------------------------
%%
dirTem='/endosome/work/bioinformatics/s171152/data/CMEmodel/St1/19';
recTem=ComRecord([dirTem filesep 'data'],[dirTem filesep 'fig'],[dirTem filesep 'init'],...
              'nRep',nPar,'deleteFiles',false);
for iPar=9:11 
    dirTem=dir(recTem.dir_all{iPar}.dir_data);
    S=load([dirTem(end).folder filesep dirTem(end).name]);
    Mpar{iPar}=S.M;
end
clear S;
%%
nStep=2000;
for iPar=1:nPar
    Mpar{iPar}.TypForce.pm.k_ModClathrin_ModFreeParticle=100;
    Mpar{iPar}.TypForce.pm.k_ModClathrin(1)=60.; %100
    Mpar{iPar}.TypForce.pm.k_ModClathrin(2)=200.; %1000
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModClathrin}.pm.n_min_adp=3;
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModClathrin}.pm.n_max_adp=3;
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModFreeParticle}.prop = {'Particle','needFollow'};
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.DlocRelax=20;
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.mu=1;
end
for iPar=1:2:nPar
    Mpar{iPar}.TypForce.pm.k_ModClathrin(1)=150.;
end

rMaxFromCtr=7.5; %5(6), 6(10), 6.5(13)
for iClathrin=19:20
dirRootSt1Tem=[dirRootSt1 filesep num2str(iClathrin)];
rec=ComRecord([dirRootSt1Tem filesep 'data'],[dirRootSt1Tem filesep 'fig'],[dirRootSt1Tem filesep 'init'],...
              'nRep',12,'deleteFiles',true);
for iPar=1:nPar          
    dirPar{iPar}=rec.dir_all{iPar};
end          
fprintf('add Clathrin...\n');
parfor iPar=1:nPar
% for iPar=1:1  
    M=Mpar{iPar};
    for iCtr=1:nCtr
    changed=false;
    while (changed==false)
        [M,changed] = add(M.mod{M.i_mod.ModClathrin},M,CtrPar{iPar}(iCtr,:),'rMaxFromCtr',rMaxFromCtr);
    end
    end
    Mpar{iPar}=M;
end

name={'ModClathrin_ModFreeParticle'};
parfor iPar=1:nPar
% for iPar=1:1  
    M=Mpar{iPar};
    idNew1=M.mod{M.i_mod.ModClathrin}.var.n_coord;
    idNew2=(M.mod{M.i_mod.ModFreeParticle}.var.n_coord-2:M.mod{M.i_mod.ModFreeParticle}.var.n_coord);
    [M] = locDyn(M.mod{M.i_mod.ModClathrin},M,name,idNew1,idNew2,'sMax',1,'update',false,'initOtherInfo',false,'addOtherInfo',true);
    Mpar{iPar}=M;
end

fprintf('relaxing. ..\n');
eStore=cell(nPar,1); %ModClathrin, ModClathrin_ModFreeParticle, ModMembrane
for iPar=1:nPar
    eStore{iPar}=[];
end
parfor iPar=1:nPar
% for iPar=1:1
    M=Mpar{iPar};
    %setup dynamics computation, indicating simulation steps, type of
    %forces included, whether adaptive time step is needed and so forth
    dyn=dynamics(nStep);
    Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle';'ModMembrane_ModSubstrate'};
%     Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle'};
    [dyn] = Preparation(dyn,M,Fname);
    i_loop_done=0;
    fig=ComPlot.PlotMod(Mpar{iPar}, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
    ComRecord.addFrame(dirPar{iPar},Mpar{iPar},fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false); 
    while (i_loop_done<nStep)
        [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',0.00000025,'plot_or_not',false,'sMax', 0.1,'iPar',iPar);
        f=M.TypForce;
        eStore{iPar}=[eStore{iPar};[f.int_V.ModClathrin,f.int_V.ModMembrane,f.int_V.ModClathrin_ModFreeParticle]];
        i_loop_done=i_loop_done+CutOff.nS;
        %plot and save results
        fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
        ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false,'iStage',2); 
        fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
    end
    Mpar{iPar}=M;
end
save([dirRootSt1Tem filesep 'matlab']);
figure;
for iPar=1:nPar
    for j=1:3
    subplot(1,3,j);
    plot(zscore(eStore{iPar}(:,j))); hold on;
    end
end
end
%%
nStep=20000;
rec=ComRecord([dirRootSt2 filesep 'data'],[dirRootSt2 filesep 'fig'],[dirRootSt2 filesep 'init'],...
              'nRep',nPar,'deleteFiles',true);
for iPar=1:nPar          
    dirPar{iPar}=rec.dir_all{iPar};
end
for iPar=1:nPar
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.k_c=100+0*iPar; %100
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.mu=100000; %500000
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.Vdh.V0=0.1; %0.2
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.f_const_std_std=0.0005;
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.DlocRelax=20;
    Mpar{iPar}.TypForce.int_stored.ModMembrane=[];
end
fprintf('relaxing. ..\n');
eStore=cell(nPar,1); %ModClathrin, ModClathrin_ModFreeParticle, ModMembrane
for iPar=1:nPar
    eStore{iPar}=[];
end
parfor iPar=1:nPar
% for iPar=10:10
    M=Mpar{iPar};
    %setup dynamics computation, indicating simulation steps, type of
    %forces included, whether adaptive time step is needed and so forth
    dyn=dynamics(nStep);
%     Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle';'ModMembrane_ModSubstrate'};
    Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle'};
    [dyn] = Preparation(dyn,M,Fname);
    i_loop_done=0;
    while (i_loop_done<nStep)
        [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',0.00000001,'plot_or_not',false,'sMax', 1);
        f=M.TypForce;
        eStore{iPar}=[eStore{iPar};[f.int_V.ModClathrin,f.int_V.ModMembrane,f.int_V.ModClathrin_ModFreeParticle]];
        i_loop_done=i_loop_done+CutOff.nS;
        %plot and save results
        fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
        ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false,'iStage',2); 
        fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
    end
    Mpar{iPar}=M;
end
figure;
for iPar=1:nPar
    for j=1:3
    subplot(1,3,j);
    plot((eStore{iPar}(:,j))); hold on;
    end
end
%%
dirMod='/home2/s171152/codes/matlab/mine/git/DanuserLab/biophysicsmodels';addpath(genpath(dirMod));
load('/endosome/work/bioinformatics/s171152/data/CMEmodel/matlab3.mat')
nStep=200000; 
k1=90; k2=100;
nameCond='WT';
% nameCond='CALM';
dirTem=[dirRootSt3 '_' nameCond '_' num2str(k1) '_' num2str(k2)];
% dirTem=dirRootSt2;
loadDone=true;
if loadDone==true
    rec=ComRecord([dirTem filesep 'data'],[dirTem filesep 'fig'],[dirTem filesep 'init'],...
              'nRep',12,'deleteFiles',false);
else
    rec=ComRecord([dirTem filesep 'data'],[dirTem filesep 'fig'],[dirTem filesep 'init'],...
              'nRep',12,'deleteFiles',true);
end
nPar=numel(rec.dir_all);
i_loop_start=cell(nPar,1);
for iPar=1:nPar
i_loop_start{iPar}=0;
end
for iPar=1:nPar
    if loadDone==true
    %----------------------------------------------------------------------
    dirInfo=dir(rec.dir_all{iPar}.dir_data);
    nIdxMax=0;
    nFile=numel(dirInfo);
    iFileRead=0;
    for iFile=1:nFile
    nChar=numel(dirInfo(iFile).name);
    for iChar=nChar:-1:1
    if strcmp(dirInfo(iFile).name(iChar),'_')
        istart=iChar;
        break;
    end
    end
    nIdx=str2double(dirInfo(iFile).name(iChar+1:end-4));
    if nIdx>nIdxMax
        nIdxMax=nIdx;
        iFileRead=iFile;
    end
    end    
    S=load([dirInfo(iFileRead).folder filesep dirInfo(iFileRead).name]);
    i_loop_start{iPar}=nIdxMax;
    Mpar{iPar}=S.M;
    end
    %----------------------------------------------------------------------
    Mpar{iPar}.TypForce.pm.k_ModClathrin_ModFreeParticle=200;
%     Mpar{iPar}.TypForce.pm.k_ModClathrin(1)=20*3; %100
%     Mpar{iPar}.TypForce.pm.k_ModClathrin(2)=200.; %1000
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.mu=50000+10000*iPar;
%     Mpar{iPar}.TypForce.otherInfo.ModClathrin_ModFreeParticle.idC_idFP(1:3:end,:)=[];
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.P=0.;
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.k_A=0.;
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.k_c=100+0*iPar; %100
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.mu=500000; %500000
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.Vdh.V0=0.1; %0.2
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.f_const_std_std=0.0005;
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.DlocRelax=20;
    Mpar{iPar}.TypForce.int_stored.ModMembrane=[];
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModFreeParticle}.prop = {'Particle'};
end
for iPar=1:nPar          
    dirPar{iPar}=rec.dir_all{iPar};
end
for iPar=1:6
    Mpar{iPar}.TypForce.pm.k_ModClathrin(1)=k1;
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.mu=500000/200*Mpar{iPar}.TypForce.pm.k_ModClathrin(1);
end
for iPar=7:12
    Mpar{iPar}.TypForce.pm.k_ModClathrin(1)=k2;
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.mu=500000/200*Mpar{iPar}.TypForce.pm.k_ModClathrin(1);
end
fprintf('relaxing. ..\n');
parfor iPar=1:nPar
% for iPar=nPar:nPar
    M=Mpar{iPar};
    %setup dynamics computation, indicating simulation steps, type of
    %forces included, whether adaptive time step is needed and so forth
    dyn=dynamics(nStep);
    Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle'};
%     Fname={'ModMembrane'};
    [dyn] = Preparation(dyn,M,Fname);
    i_loop_done=i_loop_start{iPar};
    while (i_loop_done<nStep)
        [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',0.00000005,'plot_or_not',false,'sMax', 0.5);
        i_loop_done=i_loop_done+CutOff.nS;
        %plot and save results
        fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
        ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false,'iStage',2); 
        fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
    end
    loop_done_pre{iPar}=i_loop_done;
    Mpar{iPar}=M;
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% ========================================================================
% dirTem=[dirRootSt3 '_alt'];
% dirTem=dirRootSt3;
% dirTem=[dirRootSt2 '_alt'];
% dirTem=[dirRootSt3 '_new'];
% dirTem=dirRootSt2;
dirTem=[dirRootSt3 '_WT_30_40'];
% dirTem=[dirRootSt3 '_new1_32_rev'];
% dirTem=[dirRootSt3 '_new0_5_2_5'];
% dirTem=[dirRootSt3 '_new23'];
rec=ComRecord([dirTem filesep 'data'],[dirTem filesep 'fig'],[dirTem filesep 'init'],...
              'nRep',12,'deleteFiles',false);

nPar=numel(rec.dir_all);
iFrame=cell(nPar,1);
eStore=cell(nPar,1);
pathData=cell(nPar,1);
K=cell(nPar,1);
for iPar=1:nPar
    pathData{iPar}=rec.dir_all{iPar}.dir_data;
end
parfor iPar=1:nPar
    iFrame{iPar}=[];
    eStore{iPar}=[];
    istart=0;
    %----------------------------------------------------------------------
    dirInfo=dir(pathData{iPar});
    nIdxMax=0;
    nFile=numel(dirInfo);
    iFileRead=0;
    for iFile=3:nFile
    nChar=numel(dirInfo(iFile).name);
    for iChar=nChar:-1:1
    if strcmp(dirInfo(iFile).name(iChar),'_')
        istart=iChar;
        break;
    end
    end
    nIdx=str2double(dirInfo(iFile).name(iChar+1:end-4));
    iFrame{iPar}=[iFrame{iPar};nIdx];
    S=load([dirInfo(iFile).folder filesep dirInfo(iFile).name]);
    f=S.M.TypForce;
    [f,m,~] = ModMembrane(f,S.M);
    idx=(abs(m.var.coord(:,1))<5) & (abs(m.var.coord(:,2))<5) & ((m.var.coord(:,3))<0);
%     H=m.pm.k_c*sum(0.5*(2*m.var.f.kH(idx,:)).^2.*m.var.f.A(idx,:));
    H=m.pm.k_c*mean(0.5*(2*m.var.f.kH(idx,:)).^2);
    eStore{iPar}=[eStore{iPar};[f.int_V.ModClathrin,...
                                f.int_V.ModMembrane,...
                                f.int_V.ModClathrin/S.M.mod{S.M.i_mod.ModClathrin}.var.n_coord,...
                                f.int_V.ModMembrane/S.M.mod{S.M.i_mod.ModMembrane}.var.n_coord,...
                                H,...
                                f.int_V.ModClathrin_ModFreeParticle]];
    K{iPar}=S.M.TypForce.pm.k_ModClathrin(1);
    end 
    %----------------------------------------------------------------------
end
for iPar=1:nPar
    [~,idTem]=sort(iFrame{iPar});
    eStore{iPar}=eStore{iPar}(idTem,:);
    iFrame{iPar}=iFrame{iPar}(idTem,:);
end
%%
% figure;
for iPar=1:nPar
    if iPar~=0
    subplot(1,3,1);
    plot(iFrame{iPar},eStore{iPar}(:,1)); hold on;
    subplot(1,3,2);
    plot(iFrame{iPar},eStore{iPar}(:,5),'.-'); hold on;
    subplot(1,3,3);
    plot(K{iPar},eStore{iPar}(end,1),'*'); hold on;
    end
end
%end