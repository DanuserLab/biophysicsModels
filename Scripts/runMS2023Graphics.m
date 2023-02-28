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
%%
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
%%
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
%%
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
%%
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

%%
dirRootTem='/endosome/work/bioinformatics/s171152/data/CMEmodel/St3_WT_';
nFolder=5;
dirTem=cell(nFolder,1);
for iFolder=1:nFolder
    dirTem{iFolder}=[dirRootTem num2str((iFolder*2-1)*10) '_' num2str((iFolder*2)*10)];
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
    dataToPlot{iFolder}(iPar,2)=sqrt(f.int_V.ModClathrin)/nLeg/lPL;
    %----------------------------------------------------------------------
end
end
%%
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
%% plot phase diagram histogram
TFall=[];
for iFolder=1:nFolder
    TFall=[TFall;dataToPlot{iFolder}(:,2)];
end
histogram(TFall,'normalization','pdf');xlim([0.02 0.08])
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