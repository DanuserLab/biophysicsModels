function [] = exampleMotorPolymerChain(dirRoot,dirMod,varargin)
%--------------------------------------------------------------------------
        % exampleMotorPolymerChain performs the example of the interaction
        % between a @ModBrownianMotor and @ModPolymerChain
        % input: 
        % dirRoot - root directory for data storage
        % dirMod - directory of @ModMembrane and other required objects
        % optional:
        % see variable arguments
        %   See also tetherPull, membraneFussion, exampleLamellipodia,
        %   exampleEndocytosis, redBloodCell
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/08/02
%--------------------------------------------------------------------------        
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dirRoot', @(x) ischar(x));
ip.addRequired('dirMod', @(x) ischar(x));
ip.addParameter('nStep', 2000, @isnumeric); %simulation steps for membrane relaxation
ip.addParameter('nRep', 12, @isnumeric); %repeat number
ip.addParameter('xyzLim', [-15 15;-15 15;-15 15], @isnumeric); %figure axis range
ip.addParameter('viewAng', [0 90], @isnumeric); %figure view angle
ip.parse(dirRoot,dirMod,varargin{:});
% dirRoot = '/endosome/work/bioinformatics/s171152/data/membrane/BrownianMotor';
% dirMod='/home2/s171152/codes/matlab/mine/git/DanuserLab/biophysicsmodels';addpath(genpath(dirMod));
% nStep=10000;nRep=12;xyzLim=[-15 15;-15 15;-15 15];viewAng=[0 90];
%--------------------------------------------------------------------------
nStep=ip.Results.nStep;
nRep=ip.Results.nRep;
xyzLim=ip.Results.xyzLim;
viewAng=ip.Results.viewAng;
%%
%--------------------------------------------------------------------------
rec=ComRecord([dirRoot filesep 'data'],...
              [dirRoot filesep 'fig'],...
              [dirRoot filesep 'init'],'nRep',nRep,'deleteFiles',true);
nPar=numel(rec.dir_all);

%--------------------------------------------------------------------------
%Compute MotorPolymerChain
%--------------------------------------------------------------------------
%setup parameters and variables
u=ComUnit('erg',ComUnit.nm_to_cm(10),300,ComUnit.kBT_to_erg(10,300));
pc=ModPolymerChain('dyn2','unit',u,'nChain',1,'nSubunit',20); 
p=ModFreeParticle(1,'unit',u,'coordAssigned',pc.var.coord(1,:));
b=ModBrownianMotor('flash',1,'ModFreeParticle','ModPolymerChain','dyn2','unit',u);
pc.var.coord(:,3)=pc.var.coord(:,3)-(1:pc.var.n_coord)'*0.15;

%assemble the model with the membrane object
int_info=struct('TypName',{'TypForce';'TypChemistry'},'IntList',{{{p,b},{pc,b}};{{}}});
M = model(int_info,dirMod,u,[-50 50; -50 50; -50 50],p,pc,b);
    %adjust parmeters for this particular application
    %------------------------------------------
    
    %------------------------------------------
%setup parallel variables to support parfor loop    
Mpar=cell(nPar,1);
dirPar=cell(nPar,1);
for iPar=1:nPar
    Mpar{iPar}=M;
    dirPar{iPar}=rec.dir_all{iPar};
end
%%
%--------------------------------------------------------------------------
%dynamics
fprintf('computing ...\n');
% parfor iPar=1:nPar
for iPar=1:1    
    M=Mpar{iPar};
    %setup dynamics computation, indicating simulation steps, type of
    %forces included, whether adaptive time step is needed and so forth
    dyn=dynamics(nStep);
%     Fname={'ModBrownianMotor_ModPolymerChain'};
    Fname={'ModPolymerChain'};
    [dyn] = Preparation(dyn,M,Fname);
    i_loop_done=0;
    %stop simulation at every Tcutoff to save results
    Tcutoff=0.1;
    while (i_loop_done<nStep)
        [M,CutOff] = TimeEval(dyn,M,Fname,'tCutoff',Tcutoff,'sMax',0.02,'plot_or_not', true);
        i_loop_done=i_loop_done+CutOff.nS;
        fprintf('loop finished %d at par %d\n', i_loop_done,iPar);
        %plot and save results
        fig=ComPlot.PlotMod(M, 'xyzLim', xyzLim,'viewAng',viewAng);
        ComRecord.addFrame(dirPar{iPar},M,fig,i_loop_done,'rec_fig_only',false,'update_mod_only',false);
        %--------------------------------------------------------------------------
    end
    ComPlot.PlotMod(M, 'xyzLim', xyzLim,'viewAng',viewAng);
Mpar{iPar}=M;
%makeMovie(rec,'DelayTime', 0.25);
fprintf('finishded at i_par %d \n',iPar);
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
end