%==========================================================================
% This script includes all the examples studied in the 2022 manuscript. 
% Applications are separated into sections by double %% signs 
% and can be run individually with 'ctrl'+'enter'
% see the documentation of individual application for more detail
%   See also runMS2022Graphics
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/25
%==========================================================================
%% add library
dirMod='/home2/s171152/codes/matlab/mine/git/DanuserLab/biophysicsmodels'; 
addpath(genpath(dirMod));
%% ========================================================================
% Computation
%==========================================================================
%% Red blood cell
dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/Fig3_RBC';
BiophysicsApp.redBloodCell(dirRoot,dirMod,'nRep', 12); %repeat 12 times
%% Fusion
dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/Fig3_fussion';
BiophysicsApp.membraneFussion(dirRoot,dirMod,'nRep', 12);
%% Filopodia
dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/test';
BiophysicsApp.exampleFilopodia(dirRoot,dirMod,'nRep', 12);
%% Lamellipodia
dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/Fig4_Lamellipodia';
BiophysicsApp.exampleLamellipodia(dirRoot,dirMod,'nRep', 12);
%% Endocytosis
dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/Fig4_Endocytosis';
BiophysicsApp.exampleEndocytosis(dirRoot,dirMod,'nRep', 12);
%% Tether pulling
dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/Fig5';
BiophysicsApp.tetherPull(dirRoot,dirMod,'kDiffAlt', [0.025,0.05,0.1,0.15,0.2,0.25],'reloadInit',true,'idInitNormal',7,'idInitExtraMem',10,...
                        'reloadInitExtraMem',true,'reloadComputed',true,'reloadAnalysis',true,'r_threAlt',0.1);