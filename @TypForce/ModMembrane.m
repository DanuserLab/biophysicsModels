function [f,abnormal_length,m,Vtot] = ModMembrane(f,mod,varargin)
%--------------------------------------------------------------------------
        % ModMembrane performs the computation of the forces involved in
        % @ModeMembrane, including bending force, internal force, pressure
        % force and tension force
        % input: 
        % f - @TypForce
        % mod - @model object
        % output:
        % abnormal_length - whether too long or too short edge appear
        % m - @ModeMembrane object with updated force data
        % Vtot - total potential value of m
        % optional:
        % see variable arguments
        %   See also ModMembrane_ModSubstrate
        % Author: Xinxin Wang, Danuser Lab
        % email: wangxinxin8627@gmail.com
        % date: 2022/07/19
%--------------------------------------------------------------------------
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addParameter('mex', [], @isobject);
ip.addParameter('Pflat', false, @islogical);
ip.addParameter('split', [], @isstruct);
ip.addParameter('merge', [], @isstruct);
ip.parse(f,mod,varargin{:});
%----------------------------------------------------------------------------------------
m=mod.mod{mod.i_mod.ModMembrane};
%----------------------------------------------------------------------------------------
mex_avail=true;
%----------------------------------------------------------------------------------------
if isempty(f.int_stored.ModMembrane)
    f.int_stored.ModMembrane = ModMembrane8Const(f,m);
end
%----------------------------------------------------------------------------------------
if m.pm.remeshScheme==0
    Vpm=m.pm.Vdw;
else
    Vpm=m.pm.Vdh;
end
%----------------------------------------------------------------------------------------
d = (m.var.coord(m.var.edge_all(m.var.id_on_edg,2),:) - m.var.coord(m.var.edge_all(m.var.id_on_edg,1),:));
r = sqrt(sum(d.^2,2));
u = d./r;
i_shift=f.int_stored.ModMembrane.rn(1)/m.pm.dr-1;
i = floor(r/m.pm.dr+0.5)-i_shift;

f_edg=f.int_stored.ModMembrane.fn(i).*u;
f.int_const.ModMembrane=zeros(m.var.n_coord,3);
Vtot=sum(f.int_stored.ModMembrane.Vn(i));
for i_coord = 1:numel(m.var.id_on_coord)
    f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:) = f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:)-sum(f_edg(m.var.edge_all(m.var.id_on_edg,1)==m.var.id_on_coord(i_coord),:),1);
    f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:) = f.int_const.ModMembrane(m.var.id_on_coord(i_coord),:)+sum(f_edg(m.var.edge_all(m.var.id_on_edg,2)==m.var.id_on_coord(i_coord),:),1);
end
%----------------------------------------------------------------------------------------
if m.pm.remeshScheme==0
    id_tem1 = r>Vpm.rl_max; 
    n_tem1 = length(id_tem1(id_tem1));
    if (n_tem1 == 0) 
    abnormal_length=false;
    else
    abnormal_length=true;
    end
else
id_tem1 = r<Vpm.rl_min;  id_tem2 = r>Vpm.rl_max; 
n_tem1 = length(id_tem1(id_tem1));
n_tem2 = length(id_tem2(id_tem2));
if (n_tem1 == 0) && (n_tem2 == 0)
    abnormal_length=false;
else
    abnormal_length=true;
end
end
%----------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------
f.int_rand.ModMembrane = randn(m.var.n_coord,3)*sqrt(2*m.pm.mu*m.pm.kBT*m.pm.dt)/m.pm.dt/m.pm.mu;

if mex_avail == false
%%
[m] = f.ModMembrane8Helfrich(m,'idx',m.var.id_on_coord);
coord_org = m.var.coord;
H_org = m.var.f.H;
Vtot=Vtot+H_org;
K_org=m.var.f.K;
kH_org = m.var.f.kH;
A_org = m.var.f.A;
f.int_comp.ModMembrane = zeros(m.var.n_coord,3);
m.var.f.dH = m.var.f.dH*0;
dr_H=0.00001;
for i_on = 1:m.var.n_on_coord
    i=m.var.id_on_coord(i_on);
    for k=1:3
       m.var.coord = coord_org; m.var.coord(i,k) = m.var.coord(i,k) + dr_H;
       idx = [i,m.var.j_T(i,1:m.var.val(i))];
       [~,id_tem]=intersect(idx,m.var.id_on_coord);
       idx=idx(id_tem);
       m = f.ModMembrane8Helfrich(m, 'idx',idx,'init',false);
       H_new = m.var.f.H;
%        id_tem=1:m.var.n_coord;
%        for i_tem=1:max(size(idx))
%            id_tem(id_tem==idx(i_tem))=[];
%        end
%        norm(kH_org(id_tem,:) - m.var.f.kH(id_tem,:))
%        pause
       f.int_comp.ModMembrane(i,k) = -m.pm.k_c*(H_new-H_org)/dr_H;
       m.var.f.dH(i,k) = (H_new-H_org)/dr_H;
       m.var.f.kH(idx)=kH_org(idx);
       m.var.f.A(idx)=A_org(idx);
       m.var.f.K(idx,:)=K_org(idx,:);
    end
end
%%
K_tem = m.var.f.K./m.var.f.A;
K_n=sqrt(sum(K_tem.*K_tem,2));
u_K = K_tem./K_n;

m.var.coord = coord_org;
else
    pmc=zeros(10,1);
    pmc(1) = m.pm.nt;
    pmc(2) = m.pm.dt;
    pmc(3) = m.pm.P;
    pmc(4) = m.pm.k_c;
    pmc(5) = m.pm.k_e;
    pmc(6) = m.pm.dr;
    pmc(7) = m.pm.k_V;
    pmc(8) = m.pm.k_A;
    pmc(9) = m.pm.V0;
    pmc(10) = m.pm.A0;
    pmc(11) = m.pm.nAVmean;
    pmc(12) = m.pm.k_a;
       
    j_T = m.var.j_T; j_T(isnan(j_T)) = 0;

    [f_b,f_p,u_K,kH,f_AV,V_A_Vol]=Mex.ModMembrane...
       (m.var.coord,...
        pmc,...
        m.var.edge_all,...
        m.var.face_unq,...
        j_T,...
        m.var.id_on_coord,...
        m.var.n_node',...
        m.var.T_s',...
        m.var.T_e',...
        m.var.dens);
    %%
%     fig=figure('units','normalized','outerposition',[0 0 1 1]);
%     subplot(1,2,1);
%     plot(m,'f',fig,'LineStyle','-','col',m.var.dens/max(m.var.dens));
%     subplot(1,2,2);
%     plot(m,'f',fig,'LineStyle','-','col',sum(f_AV.^2,2)/max(sum(f_AV.^2,2)));   
%     quiver3(m.var.coord(:,1),m.var.coord(:,2),m.var.coord(:,3),f_AV(:,1),f_AV(:,2),f_AV(:,3)); hold on;
%%
    m.var.f.kH=kH;
    m.var.f.u_K=u_K;
    m.var.f.f_AV=f_AV;
    if ip.Results.Pflat==true
    f_p(m.var.coord(:,3)<0.0001,:)=0;
    end
    f.int_comp.ModMembrane=f_b+f_p+f_AV;
    Vtot=Vtot+V_A_Vol(1);
    m.var.f.A=V_A_Vol(2);
    m.var.f.V=V_A_Vol(3);
end
f.int_V.ModMembrane=Vtot;
f.int_tot.ModMembrane=cell(1,1);
f.int_tot.ModMembrane{1}=f.int_const.ModMembrane+f.int_comp.ModMembrane+f.int_rand.ModMembrane;