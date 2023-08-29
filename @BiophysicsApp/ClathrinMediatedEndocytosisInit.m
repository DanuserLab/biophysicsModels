function [] = ClathrinMediatedEndocytosisInit(dirRoot,dirMod,varargin)
%--------------------------------------------------------------------------
%         ClathrinMediatedEndocytosisInit performs the initiation for
%         ClathrinMediatedEndocytosis
%         input: 
%         dirRoot - root directory for data storage
%         dirMod - directory of @ModMembrane and other required objects
%         optional:
%         see variable arguments
%           See also tetherPull, membraneFussion, exampleLamellipodia,
%           exampleEndocytosis, redBloodCell
%         Author: Xinxin Wang, Danuser Lab
%         email: wangxinxin8627@gmail.com
%         date: 2023/08/27
%--------------------------------------------------------------------------        
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dirRoot', @(x) ischar(x));
ip.addRequired('dirMod', @(x) ischar(x));
ip.addParameter('nStep', [10000 2000], @isnumeric); 
%1: simulation steps for after adding CALM control points; 
%2: simulation steps for relaxation after each clathrin added
ip.addParameter('nRep', 4, @isnumeric); %repeat number, repeat 4 times for each initial curvature: flat, shallow and high
ip.addParameter('xyzLim', [-20 20;-20 20;-20 20], @isnumeric); %mesh range
ip.addParameter('xyzLimPlot', [-10 10;-10 10;-15 15], @isnumeric); %mesh range
ip.addParameter('viewAng', [0 0], @isnumeric); %figure view angle
ip.addParameter('Tcutoff', 0.00000005, @isnumeric); %time cutoff for when to exit and plot frames
ip.addParameter('loadRes', true, @islogical); 
ip.addParameter('nCL', 20, @isnumeric);
% ip.addParameter('idCALM', [1 5 9], @isnumeric); 
ip.parse(dirRoot,dirMod,varargin{:});
%--------------------------------------------------------------------------
nRep=ip.Results.nRep;
xyzLim=ip.Results.xyzLim;
viewAng=ip.Results.viewAng;
xyzLimPlot=ip.Results.xyzLimPlot;
Tcutoff=ip.Results.Tcutoff;
loadRes=ip.Results.loadRes;
nCL=ip.Results.nCL;
%--------------------------------------------------------------------------
dirRoot=ip.Results.dirRoot;
dirMod=ip.Results.dirMod;
% dirRoot = '/endosome/work/bioinformatics/s171152/data/CMEmodel';
% dirMod='/home2/s171152/codes/matlab/mine/git/DanuserLab/biophysicsmodels';addpath(genpath(dirMod));
%--------------------------------------------------------------------------
%% initial parameters
%%
% Initialization
%--------------------------------------------------------------------------
%setup parameters and variables
% unit: 10nm, 1kBT as natural units 
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
% membrane
m=ModMembrane(true,4,0,'unit',u);
% control points for pulling membrane as init conditions
s=ModSubstrate(0.01,13);
% clathrin
c=ModClathrin('unit',u);
% later as AP2-adaptor proteins to connect clathrin and membrane
p=ModFreeParticle(0,'unit',u);
p.prop = {'Particle','needFollow'};

%introducing various interactions to be included, for example: {m,s} means ctrl points interact with membrane
int_info=struct('TypName',{'TypForce';'TypChemistry'},'IntList',{{{c,p},{m,s},{p,m}};{{}}});
% assemble the model with the objects
M = model(int_info,dirMod,u,xyzLim,m,s,c,p);
%% computation
% other related parameter setting
%--------------------------------------------------------------------------
% osmotic pressure PV
P=[1];
% membrane tension kA(A-A0)^2
kA=[8];
% setup path for saving results, iCALM indicates 3 different initial
% curvature: flat, shallow and high
dirFolder=[dirRoot filesep 'CALMinit']; mkdir(dirFolder);
if loadRes==true
    S=load([dirFolder filesep 'rec.mat']);
    rec=S.rec;
    clear S;
else
    rec=ComRecord([dirFolder filesep 'data'],[dirFolder filesep 'fig'],[dirFolder filesep 'init'],...
              'nRep',nRep,'deleteFiles',true,...
              'paraAlt',struct('P',P,'kA',kA,'iCALM',[0 1 2]));
    save([dirFolder filesep 'rec.mat'],'rec');
end
nPar=numel(rec.dir_all); %total number of workers, defult: 12
%--------------------------------------------------------------------------
% get control points' positions
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
% adjust membrane parameters, see ModMembrane for more info
%--------------------------------------------------------------------------
M.mod{M.i_mod.ModMembrane}.pm.A0=4*pi*abs(minZ)^2*0.9; % here, introduce some tension to begin
% M.mod{M.i_mod.ModMembrane}.pm.V0=4/3*pi*abs(minZ)^3*1.1;
M.mod{M.i_mod.ModMembrane}.pm.k_V=0;
M.mod{M.i_mod.ModMembrane}.pm.mu=500000; %500000
M.mod{M.i_mod.ModMembrane}.pm.Vdh.V0=0.1; %lipid reorganization
M.mod{M.i_mod.ModMembrane}.pm.f_const_std_std=0.0005;
M.mod{M.i_mod.ModMembrane}.pm.DlocRelax=20;
%--------------------------------------------------------------------------
%setup parallel variables to support parfor loop    
Mpar=cell(nPar,1);
dirPar=cell(nPar,1);
for iPar=1:nPar
    Mpar{iPar}=M;
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.P=rec.dir_all{iPar}.paraVal(1);
    Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.k_A=rec.dir_all{iPar}.paraVal(2);
    Mpar{iPar}.TypForce.pm.k_ModMembrane_ModSubstrate=10;  % for membrane-ctrl point force
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
if loadRes==true
    fprintf('CALM result reloading... %d \n',iPar);
    for iPar=1:nPar
       dirTem=dir(rec.dir_all{iPar}.dir_data);
       S=load([dirTem(end).folder filesep dirTem(end).name]);
       Mpar{iPar}=S.M;
       clear S;
    end
else
%Computation
parfor iPar=1:nPar
%     delete([dirPar{iPar}.dir_data filesep '*.*']);
%     delete([dirPar{iPar}.dir_fig filesep '*.*']);
    i_loop_done=0;
    fig=ComPlot.PlotMod(Mpar{iPar}, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
    ComRecord.addFrame(dirPar{iPar},Mpar{iPar},fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false); 
end

nStep=ip.Results.nStep(1);
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
fprintf('CALM finishded at i_par %d \n',iPar);
end
end
%==========================================================================
% adding clathrin
%==========================================================================
nPar=numel(rec.dir_all);
dirFolder=[dirRoot filesep 'CLinit']; mkdir(dirFolder);
% if loadRes==true
%     S=load([dirFolder filesep 'rec.mat']);
%     rec=S.rec;
%     clear S;
% else
    recCL=ComRecord([dirFolder filesep 'data'],[dirFolder filesep 'fig'],[dirFolder filesep 'init'],'deleteFiles',true);    
    for iPar=1:nPar
        dirPar{iPar}=recCL.dir_all{iPar};
    end
%     save([dirFolder filesep 'rec.mat'],'rec');
% end
% dirPar=cell(nPar,1);
% S=load([dirRootSt1Tem filesep 'FlatInit.mat']);
% for iPar=1:4         
%     Mpar{iPar}=S.M;
%     dirPar{iPar}=rec.dir_all{iPar};
% end
% S=load([dirRootSt1Tem filesep 'LowInit.mat']);
% for iPar=5:8
%     Mpar{iPar}=S.M;
%     dirPar{iPar}=rec.dir_all{iPar};
% end
% S=load([dirRootSt1Tem filesep 'HighInit.mat']);
% for iPar=9:12
%     Mpar{iPar}=S.M;
%     dirPar{iPar}=rec.dir_all{iPar};
% end
%--------------------------------------------------------------------------
% assign centers to recruit clathrin, currently only support one center per
% sphere of membrane
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

%%
nStep=ip.Results.nStep(2);
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

for iClathrin=2:nCL
    %%---------------------------------------------------------------------
    dirFolder=[dirRoot filesep 'CLadding' filesep num2str(iClathrin)]; mkdir(dirFolder);
    recAdd=ComRecord([dirFolder filesep 'data'],[dirFolder filesep 'fig'],[dirFolder filesep 'init'],'deleteFiles',true);    
    for iPar=1:nPar
        dirPar{iPar}=recAdd.dir_all{iPar};
    end
    %%---------------------------------------------------------------------
    if iClathrin<=6
        rMaxFromCtr=5; %5(6), 6(10), 6.5(13)
    elseif iClathrin<=10
        rMaxFromCtr=6; 
    elseif iClathrin<=13
        rMaxFromCtr=6.5; 
    else
        rMaxFromCtr=7.5; 
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
        [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',Tcutoff*5,'plot_or_not',false,'sMax', 0.1,'iPar',iPar);
        f=M.TypForce;
%         eStore{iPar}=[eStore{iPar};[f.int_V.ModClathrin,f.int_V.ModMembrane,f.int_V.ModClathrin_ModFreeParticle]];
        i_loop_done=i_loop_done+CutOff.nS;
        %plot and save results
        fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
        ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false); 
        fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
    end
    Mpar{iPar}=M;
end
% save([dirRootSt1Tem filesep 'matlab']);
% figure;
% for iPar=1:nPar
%     for j=1:3
%     subplot(1,3,j);
%     plot(zscore(eStore{iPar}(:,j))); hold on;
%     end
% end
end
% %%
% nStep=20000;
% rec=ComRecord([dirRootSt2 filesep 'data'],[dirRootSt2 filesep 'fig'],[dirRootSt2 filesep 'init'],...
%               'nRep',nPar,'deleteFiles',true);
% for iPar=1:nPar          
%     dirPar{iPar}=rec.dir_all{iPar};
% end
% for iPar=1:nPar
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.k_c=100+0*iPar; %100
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.mu=100000; %500000
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.Vdh.V0=0.1; %0.2
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.f_const_std_std=0.0005;
%     Mpar{iPar}.mod{Mpar{iPar}.i_mod.ModMembrane}.pm.DlocRelax=20;
%     Mpar{iPar}.TypForce.int_stored.ModMembrane=[];
%     Mpar{iPar}.TypForce.pm.k_ModClathrin(1)=40;
% end
% fprintf('relaxing. ..\n');
% eStore=cell(nPar,1); %ModClathrin, ModClathrin_ModFreeParticle, ModMembrane
% for iPar=1:nPar
%     eStore{iPar}=[];
% end
% parfor iPar=1:nPar
% % for iPar=10:10
%     M=Mpar{iPar};
%     %setup dynamics computation, indicating simulation steps, type of
%     %forces included, whether adaptive time step is needed and so forth
%     dyn=dynamics(nStep);
% %     Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle';'ModMembrane_ModSubstrate'};
%     Fname={'ModClathrin';'ModMembrane';'ModClathrin_ModFreeParticle'};
%     [dyn] = Preparation(dyn,M,Fname);
%     i_loop_done=0;
%     while (i_loop_done<nStep)
%         [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',0.00000001,'plot_or_not',false,'sMax', 1);
%         f=M.TypForce;
%         eStore{iPar}=[eStore{iPar};[f.int_V.ModClathrin,f.int_V.ModMembrane,f.int_V.ModClathrin_ModFreeParticle]];
%         i_loop_done=i_loop_done+CutOff.nS;
%         %plot and save results
%         fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
%         ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false,'iStage',2); 
%         fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
%     end
%     Mpar{iPar}=M;
% end
% figure;
% for iPar=1:nPar
%     for j=1:3
%     subplot(1,3,j);
%     plot((eStore{iPar}(:,j))); hold on;
%     end
% end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------