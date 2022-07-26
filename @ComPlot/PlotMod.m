function f = PlotMod(M, varargin)
%--------------------------------------------------------------------------
        % PlotMod plots all the objects in the given @Model object M
        % input: 
        % dirAll - directory structure for data and figure
        % M - a given @Model object
        % optional:
        % see variable arguments
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/06/19
%-------------------------------------------------------------------------- 
ip = inputParser;
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addParameter('f', [], @isobject);
ip.addParameter('xyzLim', [], @isnumeric);
ip.addParameter('viewAng', [0 0], @isnumeric);
ip.parse(M, varargin{:});
%==========================================================================
f = ip.Results.f;
if isempty(f)
    f=figure;
    figure(f); hold on;
else
    figure(f); hold on;
end
xyzLim=ip.Results.xyzLim;
if isempty(xyzLim)
    xyzLim=[-10 10;-5 5;-15 5];
end
viewAng=ip.Results.viewAng;

for i=1:M.n_mod
    if strcmp(M.name{i},'ModMembrane')
        plot(M.mod{M.i_mod.ModMembrane},'f',f,'LineStyle','-','facealpha',1);
    elseif strcmp(M.name{i},'ModFreeParticle')
        plot(M.mod{M.i_mod.ModFreeParticle},'f',f,'Color',[0 1 0]);
    elseif strcmp(M.name{i},'ModClathrin')
        col_tem=zeros(M.mod{M.i_mod.ModClathrin}.var.n_coord,3); col_tem(:,2)=1;col_tem(:,3)=1;
        plot(M.mod{M.i_mod.ModClathrin},'f',f,'simple',true,'col',col_tem);
    elseif strcmp(M.name{i},'ModSubstrate')
        plot(M.mod{M.i_mod.ModSubstrate},'f',f,'col',[1 0 0]);
    else
    end
end

xlim(xyzLim(1,:));ylim(xyzLim(2,:));zlim(xyzLim(3,:));
lighting gouraud;
view(viewAng);


