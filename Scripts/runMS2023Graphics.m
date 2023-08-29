%==========================================================================
% This script adjust the graphic output for the 2023 manuscript. 
% Thie script is supposed to run after finishing runMS2022Examples
% Applications are separated into sections by double %% signs 
% and can be run individually with 'ctrl'+'enter'
% see the documentation of individual application for more detail
%   See also runMS2022Examples
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/25
%==========================================================================
%% add library
dirMod='/home2/s171152/codes/matlab/mine/git/DanuserLab/biophysicsmodels'; 
addpath(genpath(dirMod));
%% ==========================================================================
% Graphics
%==========================================================================
% plot 2 clathrin triskilia
%% 1.
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
c1=ModClathrin('unit',u);
c2=ModClathrin('unit',u);
ic=1;
c2.var.a(ic,:)=c2.var.a(ic,:)+[0. 0 1];
             c2.var.ang_a.Phi(ic)=sqrt(sum((c2.var.a(ic,:)).^2,2));
             if c2.var.ang_a.Phi(ic)>pi
                 c2.var.ang_a.Phi(ic)=2*pi-c2.var.ang_a.Phi(ic);
                 c2.var.a(ic,:)=-c2.var.a(ic,:)/norm(c2.var.a(ic,:))*c2.var.ang_a.Phi(ic);
             end
             if (c2.var.ang_a.Phi(ic)>pi) || (c2.var.ang_a.Phi(ic)<0)
                 c2.var.ang_a.Phi(ic)
                 error('wrong a')
             end         
             c2.var.coord(ic,:)=c2.var.coord(ic,:)+[2 0 0]; 
         [c2.var] = setVar(c2,true); 
fig=figure;
facealpha=0.25;
plot(c1,'f',fig,'col',[0 0 1],'facealpha',facealpha);
plot(c2,'f',fig,'col',[1 0 1],'facealpha',facealpha);
xlim([-3 5]); ylim([-4 4]); zlim([-6 2]);
axis off
view([5 0])
% ComRecord.savePlot(pwd,'triskilion','figH',fig,'PaperPosition', [0 0 5 5]);
%% 1.
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
c1=ModClathrin('unit',u);
c2=ModClathrin('unit',u);
ic=1;
c2.var.a(ic,:)=c2.var.a(ic,:)+[0. 1.2 1];
             c2.var.ang_a.Phi(ic)=sqrt(sum((c2.var.a(ic,:)).^2,2));
             if c2.var.ang_a.Phi(ic)>pi
                 c2.var.ang_a.Phi(ic)=2*pi-c2.var.ang_a.Phi(ic);
                 c2.var.a(ic,:)=-c2.var.a(ic,:)/norm(c2.var.a(ic,:))*c2.var.ang_a.Phi(ic);
             end
             if (c2.var.ang_a.Phi(ic)>pi) || (c2.var.ang_a.Phi(ic)<0)
                 c2.var.ang_a.Phi(ic)
                 error('wrong a')
             end         
             c2.var.coord(ic,:)=c2.var.coord(ic,:)+[2 0.3 -1.3]; 
         [c2.var] = setVar(c2,true); 
fig=figure;
facealpha=1;
plot(c1,'f',fig,'col',[0 0 1],'facealpha',facealpha);
plot(c2,'f',fig,'col',[1 0 1],'facealpha',1);
xlim([-3.5 3.5]); ylim([-3.5 3.5]); zlim([-5.5 1.5]);
axis off
view([30 15])
%% 2.
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
c1=ModClathrin('unit',u);
c2=ModClathrin('unit',u);
ic=1;
c2.var.a(ic,:)=c2.var.a(ic,:)+[-0.2 1 1];
             c2.var.ang_a.Phi(ic)=sqrt(sum((c2.var.a(ic,:)).^2,2));
             if c2.var.ang_a.Phi(ic)>pi
                 c2.var.ang_a.Phi(ic)=2*pi-c2.var.ang_a.Phi(ic);
                 c2.var.a(ic,:)=-c2.var.a(ic,:)/norm(c2.var.a(ic,:))*c2.var.ang_a.Phi(ic);
             end
             if (c2.var.ang_a.Phi(ic)>pi) || (c2.var.ang_a.Phi(ic)<0)
                 c2.var.ang_a.Phi(ic)
                 error('wrong a')
             end         
             c2.var.coord(ic,:)=c2.var.coord(ic,:)+[2 0 -1]; 
         [c2.var] = setVar(c2,true); 
fig=figure;
facealpha=0.25;
plot(c1,'f',fig,'col',[0 0 1],'facealpha',facealpha);
plot(c2,'f',fig,'col',[1 0 1],'facealpha',facealpha);
xlim([-3 5]); ylim([-4 4]); zlim([-6 2]);
axis off
view([0 50])
% ComRecord.savePlot(pwd,'triskilion','figH',fig,'PaperPosition', [0 0 5 5]);
%% 3.
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
m=ModMembrane(true,4,0,'unit',u);
xyzLim=[-20 20;-20 20;-20 20];viewAng=[0 0];xyzLimPlot=[-13.5 13.5;-13.5 13.5;-13.5 13.5];
int_info=struct('TypName',{'TypForce';'TypChemistry'},'IntList',{{{},{}};{{}}});
M = model(int_info,dirMod,u,xyzLim,m);
facealpha=1;
fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
axis off
% xlabel('x (10 nm)');ylabel('y (10 nm)');zlabel('z (10 nm)');
% ComRecord.savePlot(pwd,'triskilion','figH',fig,'PaperPosition', [0 0 5 5]);
%% 3.
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(1,300));
m=ModMembrane(true,4,0,'unit',u);
c2=ModClathrin('unit',u);
c2.var.a(ic,:)=c2.var.a(ic,:)+[-0.2 0.5 0.2];
             c2.var.ang_a.Phi(ic)=sqrt(sum((c2.var.a(ic,:)).^2,2));
             if c2.var.ang_a.Phi(ic)>pi
                 c2.var.ang_a.Phi(ic)=2*pi-c2.var.ang_a.Phi(ic);
                 c2.var.a(ic,:)=-c2.var.a(ic,:)/norm(c2.var.a(ic,:))*c2.var.ang_a.Phi(ic);
             end
             if (c2.var.ang_a.Phi(ic)>pi) || (c2.var.ang_a.Phi(ic)<0)
                 c2.var.ang_a.Phi(ic)
                 error('wrong a')
             end         
             c2.var.coord(ic,:)=c2.var.coord(ic,:)+[0 2 18]; 
[c2.var] = setVar(c2,true);  
xyzLim=[-20 20;-20 20;-20 20];viewAng=[0 0];xyzLimPlot=[-15 15;-15 15;-15 15];
int_info=struct('TypName',{'TypForce';'TypChemistry'},'IntList',{{{},{}};{{}}});
M = model(int_info,dirMod,u,xyzLim,m);
facealpha=1;
fig=ComPlot.PlotMod(M, 'xyzLim', xyzLimPlot,'viewAng',viewAng);
plot(c2,'f',fig,'col',[0 0 1],'facealpha',facealpha);
axis off
% xlabel('x (10 nm)');ylabel('y (10 nm)');zlabel('z (10 nm)');
% ComRecord.savePlot(pwd,'triskilion','figH',fig,'PaperPosition', [0 0 5 5]);
%% plot initial membrane morphologies: WT, CALM-, Epsin1-
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St1/1/data/rep04/mod_Stg_001Fr_00000');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',M.i_mod.ModMembrane);
xlim([-8 8]);ylim([-6 10]);zlim([-12 -4]);
axis off
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St1/1/data/rep11/mod_Stg_001Fr_00000');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',M.i_mod.ModMembrane);
xlim([-8 8]);ylim([-6 10]);zlim([-12 -4]);
axis off
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St1/1/data/rep05/mod_Stg_001Fr_00000');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',M.i_mod.ModMembrane);
xlim([-8 8]);ylim([-6 10]);zlim([-12 -4]);
axis off
%% plot initial membrane morphologies showing control points: WT, CALM-, Epsin1-
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/initControlPoints.mat');
Mpar=S.Mpar;
for iPar=1:4:12
    M=Mpar{iPar};
    ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModSubstrate]);
    xlim([-8 8]);ylim([-4 12]);zlim([-14 -6]); axis off; view([0 10]);
end
%% plot initial clathrin-membrane complex: WT, CALM-, Epsin1-
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_CALM_10_20/data/rep01/mod_Stg_002Fr_00923');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim([-8 8]);ylim([-8 8]);zlim([-14 -4]);
view([0 90]);
axis off
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_10_20/data/rep01/mod_Stg_002Fr_00018');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim([-8 8]);ylim([-8 8]);zlim([-14 -4]);
view([0 90]);
axis off
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_epsin_10_20/data/rep01/mod_Stg_002Fr_00754');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim([-8 8]);ylim([-8 8]);zlim([-14 -4]);
view([0 90]);
axis off
%% plot complex relax: WT
xyzLimPlot=[-8 8;-8 8;-18 -2];
viewAng=[0 20];
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_50_60/data/rep02/mod_Stg_002Fr_00927');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_50_60/data/rep02/mod_Stg_002Fr_78129');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_50_60/data/rep02/mod_Stg_002Fr_200534');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off
%% plot complex relax: CALM-
xyzLimPlot=[-8 8;-8 8;-20 -4];
viewAng=[0 20];
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_CALM_50_60/data/rep04/mod_Stg_002Fr_00937');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_CALM_50_60/data/rep04/mod_Stg_002Fr_78396');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_CALM_50_60/data/rep04/mod_Stg_002Fr_196132');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off
%% plot complex relax: Epsin1-
xyzLimPlot=[-8 8;-8 8;-20 -4];
viewAng=[0 20];
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_epsin_50_60/data/rep07/mod_Stg_002Fr_00719');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_epsin_50_60/data/rep07/mod_Stg_002Fr_75777');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_epsin_50_60/data/rep07/mod_Stg_002Fr_191895');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

%% phase transition path
dirRootTem='/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_'; nFolder=5;
% dirRootTem='/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_CALM_'; nFolder=3;
% dirRootTem='/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_epsin_'; nFolder=3;
dirTem=cell(nFolder,1);
for iFolder=1:nFolder
    dirTem{iFolder}=[dirRootTem num2str((iFolder*2-1)*10) '_' num2str((iFolder*2)*10)];
end

% dirTem=cell(1,1); dirTem{1}='/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_Dyn2_kD=100_1_k1=_20_30_40_old';nFolder=1;
%% phase transition data
dataToPlot=cell(nFolder,1);
for iFolder=1:nFolder
    rec=ComRecord([dirTem{iFolder} filesep 'data'],[dirTem{iFolder} filesep 'fig'],[dirTem{iFolder} filesep 'init'],...
              'nRep',12,'deleteFiles',false);
    nPar=numel(rec.dir_all);

pathData=cell(nPar,1);
K=cell(nPar,1);
for iPar=1:nPar
    pathData{iPar}=rec.dir_all{iPar}.dir_data;
end
dataToPlot{iFolder}=zeros(nPar,3);
for iPar=1:nPar
    istart=0;
    %----------------------------------------------------------------------
    dirInfo=dir(pathData{iPar});
    nIdxMax=0;
    nFile=numel(dirInfo);
    iFileRead=0;
    iFileToRead=0;
    for iFile=3:nFile
    nChar=numel(dirInfo(iFile).name);
    for iChar=nChar:-1:1
    if strcmp(dirInfo(iFile).name(iChar),'_')
        istart=iChar;
        break;
    end
    end
    nIdx=str2double(dirInfo(iFile).name(iChar+1:end-4));
    if nIdxMax<nIdx
        nIdxMax=nIdx;
        iFileToRead=iFile;
    end
    end    
    S=load([dirInfo(iFileToRead).folder filesep dirInfo(iFileToRead).name]);
    f=S.M.TypForce;
    dataToPlot{iFolder}(iPar,1)=S.M.TypForce.pm.k_ModClathrin(1); 
    dataToPlot{iFolder}(iPar,3)=iFileToRead;
    S.M.TypForce.pm.k_ModClathrin=[1 0];
    [f,V_tot,~] = ModClathrin(f,S.M);
    nLeg=0;
    for i=1:S.M.mod{S.M.i_mod.ModClathrin}.var.n_coord
       nLeg=nLeg+sum(sum(S.M.mod{S.M.i_mod.ModClathrin}.var.connect(:,1,i)>0));
    end
    nLeg=nLeg*0.5;
    x=S.M.mod{S.M.i_mod.ModClathrin}.var.coord_org;
    lPL=norm(x(52,1:3)-x(8,1:3));
    dataToPlot{iFolder}(iPar,2)=1./(sqrt(f.int_V.ModClathrin)/nLeg/lPL);
    %----------------------------------------------------------------------
end
end
%% phase transition
fig=figure;
Xscale=4;
for iFolder=1:nFolder
    colP=[0 0 1];
    x1=dataToPlot{iFolder}(1,1)/Xscale;
    x2=dataToPlot{iFolder}(7,1)/Xscale;
    y1=dataToPlot{iFolder}(1:6,2);
    y2=dataToPlot{iFolder}(7:12,2);
    H=ComPlot.notBoxPlot(y1,x1,'jitter',2); hold on;
%     set([H.mu],'color',colTem2);
    set(H.sdPtch,'FaceColor',colP,'EdgeColor','none');
    set(H.semPtch,'FaceColor',colP*0.3,'EdgeColor','none');
    H=ComPlot.notBoxPlot(y2,x2,'jitter',2); hold on; 
    set(H.sdPtch,'FaceColor',colP,'EdgeColor','none');
    set(H.semPtch,'FaceColor',colP*0.3,'EdgeColor','none');
end
xticks((10:10:100)/Xscale);
   nTick=numel(xticks);
   ticklabels={};
   for i=1:nTick
       ticklabels=[ticklabels;{num2str(i*10)}];
   end
xticklabels(ticklabels);
xlabel('k_1');
ylabel('TF');
xlim([0 110/Xscale])
%% phase transition for 3 k-values 
% fig=figure;
fig=gcf; hold on;
Xscale=4;
for iFolder=1:nFolder
    colP=[0 0 1];
    x1=dataToPlot{iFolder}(1,1)/Xscale;
    x2=dataToPlot{iFolder}(5,1)/Xscale;
    x3=dataToPlot{iFolder}(9,1)/Xscale;
    y1=dataToPlot{iFolder}(1:4,2);
    y2=dataToPlot{iFolder}(5:8,2);
    y3=dataToPlot{iFolder}(9:12,2);
    
    H=ComPlot.notBoxPlot(y1,x1,'jitter',2); hold on; 
    set(H.sdPtch,'FaceColor',colP,'EdgeColor','none');
    set(H.semPtch,'FaceColor',colP*0.3,'EdgeColor','none');
    set([H.data],'MarkerFaceColor',[0,0,1],'markerEdgeColor',[0,0,1],'MarkerSize',8,'Marker','*');

    H=ComPlot.notBoxPlot(y2,x2,'jitter',2); hold on; 
    set(H.sdPtch,'FaceColor',colP,'EdgeColor','none');
    set(H.semPtch,'FaceColor',colP*0.3,'EdgeColor','none');
    set([H.data],'MarkerFaceColor',[1,0,1],'markerEdgeColor',[1,0,1],'MarkerSize',8,'Marker','^');
    
    H=ComPlot.notBoxPlot(y3,x3,'jitter',2); hold on;
    set(H.sdPtch,'FaceColor',colP,'EdgeColor','none');
    set(H.semPtch,'FaceColor',colP*0.3,'EdgeColor','none');
    set([H.data],'MarkerFaceColor',[1,0,0],'markerEdgeColor',[1,0,0],'MarkerSize',8,'Marker','pentagram');
end
xticks((10:10:100)/Xscale);
   nTick=numel(xticks);
   ticklabels={};
   for i=1:nTick
       ticklabels=[ticklabels;{num2str(i*10)}];
   end
xlim([5/Xscale 75/Xscale]);   
xticklabels(ticklabels);
%% plot k3 sensitivity, k4, inserts
% dirTem='/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_K3=300_k1=_20_30_40/data'; xyzLimPlot=[-8 8;-8 8;-24 -12];
% dirTem='/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_K3=100_k1=_20_30_40/data'; xyzLimPlot=[-8 8;-8 8;-22 -10];
% dirTem='/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_Dyn2_kD=100_1_k1=_20_30_40_old/data'; xyzLimPlot=[-8 8;-8 8;-16 -4];
% dirTem='/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_Dyn2_kD=20_1_k1=_10_20_30/data';xyzLimPlot=[-8 8;-8 8;-19 -7];
dirTem='/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_Dyn2_kD=10_1_k1=_10_20_30/data';xyzLimPlot=[-8 8;-8 8;-19 -7];
dirInfo=dir(dirTem);
 
viewAng=[0 30];
for i=4:15
    dirTem=[dirInfo(i).folder filesep dirInfo(i).name];
    dirInfoTem=dir(dirTem);
    nFile=size(dirInfoTem,1);
    maxDate=0; iPicked=0;
    for iF=3:nFile
        if dirInfoTem(iF).datenum>maxDate
            maxDate=dirInfoTem(iF).datenum;
            iPicked=iF;
        end
    end
    S=load([dirInfoTem(iPicked).folder filesep dirInfoTem(iPicked).name]); M=S.M;
    ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
    xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off;
end
%% bimodality at k1=30
dirRootTem='/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_';
nFolder=2;
dirTem=cell(nFolder,1);
for iFolder=1:nFolder
    dirTem{iFolder}=[dirRootTem 'add' num2str(iFolder) '_30_30' ];
end
dataToPlot=cell(nFolder,1);
for iFolder=1:nFolder
    rec=ComRecord([dirTem{iFolder} filesep 'data'],[dirTem{iFolder} filesep 'fig'],[dirTem{iFolder} filesep 'init'],...
              'nRep',12,'deleteFiles',false);
    nPar=numel(rec.dir_all);

pathData=cell(nPar,1);
K=cell(nPar,1);
for iPar=1:nPar
    pathData{iPar}=rec.dir_all{iPar}.dir_data;
end
dataToPlot{iFolder}=zeros(nPar,3);
for iPar=1:nPar
    istart=0;
    %----------------------------------------------------------------------
    dirInfo=dir(pathData{iPar});
    nIdxMax=0;
    nFile=numel(dirInfo);
    iFileRead=0;
    iFileToRead=0;
    for iFile=3:nFile
    nChar=numel(dirInfo(iFile).name);
    for iChar=nChar:-1:1
    if strcmp(dirInfo(iFile).name(iChar),'_')
        istart=iChar;
        break;
    end
    end
    nIdx=str2double(dirInfo(iFile).name(iChar+1:end-4));
    if nIdxMax<nIdx
        nIdxMax=nIdx;
        iFileToRead=iFile;
    end
    end    
    S=load([dirInfo(iFileToRead).folder filesep dirInfo(iFileToRead).name]);
    f=S.M.TypForce;
    dataToPlot{iFolder}(iPar,1)=S.M.TypForce.pm.k_ModClathrin(1);
    dataToPlot{iFolder}(iPar,3)=iFileToRead;
    S.M.TypForce.pm.k_ModClathrin=[1 0];
    [f,V_tot,~] = ModClathrin(f,S.M);
    nLeg=0;
    for i=1:S.M.mod{S.M.i_mod.ModClathrin}.var.n_coord
       nLeg=nLeg+sum(sum(S.M.mod{S.M.i_mod.ModClathrin}.var.connect(:,1,i)>0));
    end
    nLeg=nLeg*0.5;
    x=S.M.mod{S.M.i_mod.ModClathrin}.var.coord_org;
    lPL=norm(x(52,1:3)-x(8,1:3));
    dataToPlot{iFolder}(iPar,2)=1./(sqrt(f.int_V.ModClathrin)/nLeg/lPL);
    %----------------------------------------------------------------------
end
end
%% bimodality
y=[];
for iFolder=1:nFolder
    colP=[0 0 1];
    y=[y;dataToPlot{iFolder}(1:12,2)];
end
histogram(y,'binwidth',2.5);
xlabel('');
ylabel('Repeats');
%% plot phase diagram insert
xyzLimPlot=[-7 7;-7 7;-19 -3];
viewAng=[0 20];
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_10_20/data/rep08/mod_Stg_002Fr_171457');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_30_40/data/rep02/mod_Stg_002Fr_200319');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_30_40/data/rep04/mod_Stg_002Fr_200893');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_30_40/data/rep10/mod_Stg_002Fr_200358');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off
%% plot k1 low
xyzLimPlot=[-8 8;-8 8;-20 -4];
viewAng=[0 20];
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_10_20/data/rep04/mod_Stg_002Fr_00018');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_10_20/data/rep04/mod_Stg_002Fr_77404');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_10_20/data/rep04/mod_Stg_002Fr_178732');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off
%% plot k1 high
xyzLimPlot=[-8 8;-8 8;-20 -4];
viewAng=[0 20];
S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_90_100/data/rep09/mod_Stg_002Fr_00839');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_90_100/data/rep09/mod_Stg_002Fr_75543');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off

S=load('/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_90_100/data/rep09/mod_Stg_002Fr_185361');
M=S.M;
ComPlot.PlotMod(M,'idModPlot',[M.i_mod.ModMembrane,M.i_mod.ModFreeParticle,M.i_mod.ModClathrin]);
xlim(xyzLimPlot(1,:));ylim(xyzLimPlot(2,:));zlim(xyzLimPlot(3,:)); view(viewAng);axis off
%% prepare for high-res animation
pathData={'/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_10_20/data/rep04';...
          '/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_90_100/data/rep09'};
nPar=numel(pathData);
pathRenamed=cell(nPar,1);
for iPar=1:nPar
    pathRenamed{iPar}=[pathData{iPar} filesep 'Renamed'];
    mkdir(pathRenamed{iPar});
end

for iPar=1:nPar
    dirInfo=dir(pathData{iPar});
    nIdxMax=0;
    nFile=numel(dirInfo);
    iFileRead=zeros(nFile,1);
    iFileToRead=0;
    %----------------------------------------------------------------------
    for iFile=4:nFile
    nChar=numel(dirInfo(iFile).name);
    for iChar=nChar:-1:1
    if strcmp(dirInfo(iFile).name(iChar),'_')
        istart=iChar;
        break;
    end
    end
    nIdx=str2double(dirInfo(iFile).name(iChar+1:end-4));
    iFileRead(iFile)=nIdx;
    if nIdxMax<nIdx
        nIdxMax=nIdx;
        iFileToRead=iFile;
    end
    end  
    %----------------------------------------------------------------------
    nNameMax=max(iFileRead)*10;
    nNameMax=numel(num2str(nNameMax));
    for iFile=4:nFile
        paraStageFormat=['%0.' num2str(nNameMax) 'd'];
        nameNew=num2str(iFileRead(iFile),paraStageFormat);
        copyfile([dirInfo(iFile).folder filesep dirInfo(iFile).name],[pathRenamed{iPar} filesep nameNew '.mat']);
    end
    %----------------------------------------------------------------------
end
%% ======================================================================== sensitivity
dir_data='/endosome/work/bioinformatics/s171152/data/CMEmodel/memSensitive';
mkdir(dir_data);
n_folder=4;
nPar=12;
for i_folder=1:n_folder
    disp(i_folder);
    iParStartEnd=[1 nPar];
    if i_folder==1
        dirTem=['/endosome/work/bioinformatics/s171152/data/CMEmodel/St3' '_WT_10_20'];
        iParStartEnd=[7 nPar];
    elseif i_folder==2
        dirTem=['/endosome/work/bioinformatics/s171152/data/CMEmodel/St3' '_WT_30_40'];
    elseif i_folder==3
        dirTem=['/endosome/work/bioinformatics/s171152/data/CMEmodel/St3' '_WT_50_60']; 
    elseif i_folder==4
        dirTem=['/endosome/work/bioinformatics/s171152/data/CMEmodel/St3' '_WT_70_80'];
    end

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
parfor iPar=iParStartEnd(1):iParStartEnd(2)
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
    k_ModClathrin_save=S.M.TypForce.pm.k_ModClathrin;
%     [f,m,~] = ModMembrane(f,S.M);
    m=S.M.mod{S.M.i_mod.ModMembrane};
    [m] = ModMembrane8Helfrich(f,m,'init',true);
    idx=(abs(m.var.coord(:,1))<5) & (abs(m.var.coord(:,2))<5) & ((m.var.coord(:,3))<0);
%     H=m.pm.k_c*sum(0.5*(2*m.var.f.kH(idx,:)).^2.*m.var.f.A(idx,:));
    H=m.pm.k_c*mean(0.5*(2*m.var.f.kH(idx,:)).^2.*m.var.f.A(idx));
    Htot=m.pm.k_c*sum(0.5*(2*m.var.f.kH(idx,:)).^2.*m.var.f.A(idx));
    
    f=S.M.TypForce;
    S.M.TypForce.pm.k_ModClathrin=[1 0];
    [f,V_tot,~] = ModClathrin(f,S.M);
    nLeg=0;
    for i=1:S.M.mod{S.M.i_mod.ModClathrin}.var.n_coord
       nLeg=nLeg+sum(sum(S.M.mod{S.M.i_mod.ModClathrin}.var.connect(:,1,i)>0));
    end
    nLeg=nLeg*0.5;
    x=S.M.mod{S.M.i_mod.ModClathrin}.var.coord_org;
    lPL=norm(x(52,1:3)-x(8,1:3));
    TF=1/(sqrt(f.int_V.ModClathrin)/nLeg/lPL);
    [f,m,Vtot] = ModMembrane(f,S.M);
    fmem=f.int_comp.ModMembrane;
%     fmem=f.int_tot.ModMembrane{1};
    
    eStore{iPar}=[eStore{iPar};[f.int_V.ModClathrin,...
                                f.int_V.ModMembrane,...
                                f.int_V.ModClathrin/S.M.mod{S.M.i_mod.ModClathrin}.var.n_coord,...
                                f.int_V.ModMembrane/S.M.mod{S.M.i_mod.ModMembrane}.var.n_coord,...
                                H,...
                                Htot,...
                                f.int_V.ModClathrin_ModFreeParticle,...
                                max(vecnorm(fmem(idx,:),2,2)),...
                                mean(vecnorm(fmem(idx,:),2,2)),...
                                std(vecnorm(fmem(idx,:),2,2)),...
                                TF,...
                                k_ModClathrin_save(1)]];
    K{iPar}=S.M.TypForce.pm.k_ModClathrin(1);
    end 
    %----------------------------------------------------------------------
end
for iPar=iParStartEnd(1):iParStartEnd(2)
    [~,idTem]=sort(iFrame{iPar});
    eStore{iPar}=eStore{iPar}(idTem,:);
    iFrame{iPar}=iFrame{iPar}(idTem,:);
end
   save([dir_data filesep 'iParStartEnd' num2str(i_folder)],'iParStartEnd');
   save([dir_data filesep 'iFrame' num2str(i_folder)],'iFrame');
   save([dir_data filesep 'eStore' num2str(i_folder)],'eStore');
end
%%
figure;
iStart=1;
iEnd=150;
nFigX=1;
nFigY=2;
dataAll=[];
for i_folder=1:n_folder
    load([dir_data filesep 'iParStartEnd' num2str(i_folder)]);
    load([dir_data filesep 'iFrame' num2str(i_folder)]);
    load([dir_data filesep 'eStore' num2str(i_folder)]);
%     subplot(nFigX,nFigY,1);
for iPar=iParStartEnd(1):iParStartEnd(2)
%     scatter(eStore{iPar}(iStart:iEnd,6),eStore{iPar}(iStart:iEnd,11),'r','filled','SizeData',20); alpha(.1);hold on;
%     plot(iFrame{iPar}(iStart:iEnd),eStore{iPar}(iStart:iEnd,6));hold on;
    dataAll=[dataAll;[eStore{iPar}(iStart:iEnd,6),eStore{iPar}(iStart:iEnd,11),eStore{iPar}(iStart:iEnd,12)]];
end
end
dotCol=zeros(size(dataAll,1),3);
idAll=(1:size(dataAll,1))';
idx=dataAll(:,3)==20;
dotCol(idAll(idx),1)=1;
idx=dataAll(:,3)==30;
dotCol(idAll(idx),2)=1;
idx=dataAll(:,3)==40;
dotCol(idAll(idx),3)=1;
idx=dataAll(:,3)==50;
dotCol(idAll(idx),1:2)=1;
idx=dataAll(:,3)==60;
dotCol(idAll(idx),2:3)=1;
idx=dataAll(:,3)==70;
dotCol(idAll(idx),1)=1; dotCol(idAll(idx),3)=1;
scatter(dataAll(:,1),dataAll(:,2),20,dotCol,'filled'); alpha(.2);hold on;
%%
scatter(dataAll(:,1),dataAll(:,2),20,dotCol,'filled'); alpha(.2);hold on;
hist=[];
hist_coord=[];
dT=2.5;

for T=15:dT:35
    hist_coord=[hist_coord; [T T+dT]];
    idx=(dataAll(:,2)>=T) & (dataAll(:,2)<T+dT);
    hist=[hist; numel(idx(idx==true))];
end

nBS=100;
nSamp=floor(min(hist)*0.5);
dataBS=[];
for iBS=1:nBS
    for T=15:dT:35
    idx=(dataAll(:,2)>=T) & (dataAll(:,2)<T+dT);
    Hin=dataAll(idx,1);
    Hbs=randsample(Hin,nSamp);
    dataBS=[dataBS;[median(Hbs) (2*T+dT)*0.5]];
    plot(median(Hbs), (2*T+dT)*0.5,'.'); hold on;
    end
end
%%
dataBS_reduced=[dataAll(:,1),dataAll(:,2)];
% idx=(dataBS_reduced(:,1)>35);
% dataBS_reduced(idx,:)=[];
% dataBS_reduced=dataBS;
% dataBS_reduced(nBS*4+1:nBS*5,:)=[];
xShift=-min(dataBS_reduced(:,2));
xScale=0.1;
x=(dataBS_reduced(:,2)+xShift)*xScale;
yScale=-0.001;
yShift=-3.5;
y=(dataBS_reduced(:,1))*yScale+yShift;
% x=dataBS_reduced(:,1);
% y=dataBS_reduced(:,2);
% f = fit(x,y,'exp1');
% f = fit(x,y,'poly4');
% f = fit(x,y,'gauss1');
% gaussEqn = 'a*exp(-((x-b)/c)^2)+d'; startPoints = [1.5 2 1 -3];
% f = fit(x,y,gaussEqn,'Start', startPoints);
expEqn = 'a*exp(b*x^1)+c'; startPoints = [2 -1.7 0.5];
f = fit(x,y,expEqn,'Start', startPoints);

ci=confint(f);

x_srt=sort(x);
yFit=f.a*exp(f.b*x_srt)+f.c;
y1=ci(1,1)*exp(ci(1,2)*x_srt)+ci(1,3);
y2=ci(2,1)*exp(ci(2,2)*x_srt)+ci(2,3);
% plot(x_srt,yFit,'-'); hold on;
% plot(x_srt,y1,'--'); hold on;
% plot(x_srt,y2,'--'); hold on;

x2 = [x_srt', fliplr(x_srt')];
x2=x2/xScale-xShift;
x=x/xScale-xShift;
x_srt=x_srt/xScale-xShift;
inBetween = [y2', fliplr(y1')];
inBetween=(inBetween-yShift)/yScale;
yFit=(yFit-yShift)/yScale;
y=(y-yShift)/yScale;
% scatter(y,x,'r','filled','SizeData',5); alpha(.1);hold on; 
scatter(flipud(dataAll(:,1)),flipud(dataAll(:,2)),2,dotCol,'filled'); alpha(.3);hold on;
p=patch(inBetween,x2, 'g');
p.FaceColor=[0 1 1];
p.FaceAlpha=0.9;
p.LineStyle='none';
plot(yFit,x_srt,'--','color',[0 0 1],'linewidth',1.); hold on;
xlim([500 5000]);
xlabel('E_{mem} ({\kappa})');
xticks([500 1500 2500 3500 4500]);
xticklabels({'5','15','25','35','45'});
% x_srt=sort(x);
% p22 = predint(f,x_srt,0.95,'observation','off');hold on;
% plot(x_srt,p22,'m--');
%%
x=(dataBS(:,2)-min(dataBS(:,2)))*0.1;
y=-(dataBS(:,1))*0.001;
% x=(dataBS_reduced(:,1)-min(dataBS_reduced(:,1)))*0.1;
% y=(dataBS_reduced(:,2))*0.001;
plot(x,y,'.'); hold on;
tbl = table(x, y);
modelfun = @(b,x) b(1) + b(2) * exp(b(3)*x(:, 1));  
beta0 = [-1, -1, -0.1]; 
opts = statset('Display','iter','MaxIter',10000,'TolFun',1e-20);
mdl = fitnlm(tbl, modelfun, beta0,'Options',opts);
coefficients = mdl.Coefficients{:, 'Estimate'};
% Create smoothed/regressed data using the model:
yFitted = coefficients(1) + coefficients(2) * exp(coefficients(3)*x);
plot(x, yFitted, 'r*');