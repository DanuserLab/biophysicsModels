function [mod,changed] = ModClathrin8init(ch,mod,r_init,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('ch', @(x) isa(x,'TypChemistry'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addRequired('r_init', @(x) isnumeric(x));
ip.addParameter('dr_adp', 1, @isnumeric);
ip.addParameter('update', false, @islogical);
ip.parse(ch,mod,r_init,varargin{:});
%----------------------------------------------------------------------------------------
i_mod=[mod.i_mod.ModClathrin,mod.i_mod.ModMembrane,mod.i_mod.ModMemAdapter];  %1: clathrin, 2: membrane, 3: membrane adapter
dr_adp=ip.Results.dr_adp;
update=ip.Results.update;
%----------------------------------------------------------------------------------------
%%
changed=false;
id_feet=[49 50 51];

shift_xyz=rand(100,3)*norm(mod.mod{1}.var.coord_org(1,1:3)-mod.mod{1}.var.coord_org(2,1:3))*5;
% shift_xyz=zeros(1,3);

n_xyz=size(shift_xyz,1);
d_all=inf(mod.mod{i_mod(1)}.pm.n_test_all,n_xyz);
n_on_mem=zeros(mod.mod{i_mod(1)}.pm.n_test_all,n_xyz);
i_on_mem=false(mod.mod{i_mod(1)}.pm.n_test_all,3,n_xyz);
i_c_mem=zeros(mod.mod{i_mod(1)}.pm.n_test_all,3,n_xyz);
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
for i_xyz=1:n_xyz
%--------------------------------------------------------------------------
    coord_new=r_init+shift_xyz(i_xyz,:);
%--------------------------------------------------------------------------
%%
for i_test=1:mod.mod{i_mod(1)}.pm.n_test_all
   %-----------------------------------------------------------------------
   foot=mod.mod{i_mod(1)}.get_r_from_a(mod.mod{i_mod(1)}.var.coord_org(id_feet,1:3),3,mod.mod{i_mod(1)}.pm.O_test_all(:,:,i_test))+coord_new;
   %-----------------------------------------------------------------------
   %occupied=mod.mod{i_mod(3)}.var.occupied;
   d_all(i_test,i_xyz)=0;
   n_adp_bound=0;
   for i_foot=1:3
       d=sqrt(sum((foot(i_foot,:)-mod.mod{i_mod(3)}.var.coord).^2,2));
       %d(occupied)=dr_adp+1;
       [d_min,id_min]=min(d);
       if d_min<dr_adp
       n_on_mem(i_test,i_xyz)=n_on_mem(i_test,i_xyz)+1;
       i_c_mem(i_test,i_foot,i_xyz)=id_min;
       i_on_mem(i_test,i_foot,i_xyz)=true;
       d_all(i_test,i_xyz)=d_all(i_test,i_xyz)+d_min;
       %occupied(id_min)=true;
       n_adp_bound=n_adp_bound+1;
       end
   end
   if n_adp_bound~=mod.mod{i_mod(1)}.pm.n_init_adp
       d_all(i_test,i_xyz)=inf;
   end
   %-----------------------------------------------------------------------
%    [id_reps_m_leg] = ModClathrin_ModMembrane8getIDrepSingle(mod.TypForce,mod,coord_new,mod.mod{i_mod(1)}.pm.O_test_all(:,:,i_test));
%    if isempty(id_reps_m_leg)
%        allow_add=true;
%    else
% %        id_tem = sum(abs(id_reps_m_c_leg(:,2)-i_c),2)==0;
% %        if (isempty(n_in_mem(id_tem))) || (sum(n_in_mem(id_tem))<3)
% %            allow_add=true;
% %        else
% %            allow_add=false;
% %        end
%        allow_add=false;
%    end
   
%    if  (allow_add==false) %(n_on_mem(i_test,i_xyz)>=1)
%        d_all(i_test,i_xyz)=inf;
%    end
   
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------    
%%
if numel(d_all(~isinf(d_all)))>0
[~,id_xyz]=min(d_all,[],1);
id_keep=zeros(2,1);
min_tem=inf;
for i_xyz=1:n_xyz
    if min_tem>d_all(id_xyz(i_xyz),i_xyz)
       min_tem=d_all(id_xyz(i_xyz),i_xyz);
       id_keep=[id_xyz(i_xyz),i_xyz];
    end
end
%--------------------------------------------------------------------------
coord_new=r_init+shift_xyz(id_keep(2),:);
if update==false
    mod.mod{i_mod(1)}.var.n_coord=1;
    mod.mod{i_mod(1)}.var.connect=zeros(3,2,1);
    mod.mod{i_mod(1)}.var.coord=coord_new;
    mod.mod{i_mod(1)}.var.a=mod.mod{i_mod(1)}.pm.a_test_all(id_keep(1),:);
    mod.mod{i_mod(1)}.var.ang_a.Phi=mod.mod{i_mod(1)}.pm.ang_a_test_all(id_keep(1),1);
else
    mod.mod{i_mod(1)}.var.n_coord=mod.mod{i_mod(1)}.var.n_coord+1;
    mod.mod{i_mod(1)}.var.connect=cat(3,mod.mod{i_mod(1)}.var.connect,zeros(3,2));
    mod.mod{i_mod(1)}.var.coord=[mod.mod{i_mod(1)}.var.coord;coord_new];
    mod.mod{i_mod(1)}.var.a=[mod.mod{i_mod(1)}.var.a;mod.mod{i_mod(1)}.pm.a_test_all(id_keep(1),:)];
    mod.mod{i_mod(1)}.var.ang_a.Phi=[mod.mod{i_mod(1)}.var.ang_a.Phi;mod.mod{i_mod(1)}.pm.ang_a_test_all(id_keep(1),1)];
end

mod.mod{i_mod(1)}.var.O=mod.mod{i_mod(1)}.Omega(mod.mod{i_mod(1)}.var.a,mod.mod{i_mod(1)}.var.ang_a.Phi);
mod.mod{i_mod(1)}.var.E=mod.mod{i_mod(1)}.getE(mod.mod{i_mod(1)}.var.a,mod.mod{i_mod(1)}.var.ang_a.Phi);
mod.mod{i_mod(1)}.var.bound(mod.mod{i_mod(1)}.var.n_coord,:)=i_on_mem(id_keep(1),:,id_keep(2));

for i_leg_tem=1:3
    if i_c_mem(id_keep(1),i_leg_tem,id_keep(2)) > 0
      mod.mod{i_mod(3)}.var.id_ModClathrin(i_c_mem(id_keep(1),i_leg_tem,id_keep(2)),:)=[mod.mod{i_mod(1)}.var.n_coord i_leg_tem];
      mod.mod{i_mod(3)}.var.occupied(i_c_mem(id_keep(1),i_leg_tem,id_keep(2)))=true;
    end
end
changed=true;
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


% fig=figure;
% plot(mod.mod{i_mod},'col',[],'f',fig);
% hold on;
% r_tem1=mod.mod{i_mod}.get_r_from_a([0,0,1],1,mod.mod{i_mod}.var.O(:,:,i_c));
% quiver3(mod.mod{i_mod}.var.coord(i_c,1),mod.mod{i_mod}.var.coord(i_c,2),mod.mod{i_mod}.var.coord(i_c,3),r_tem1(1),r_tem1(2),r_tem1(3),...
%         10,'linewidth',2);
% hold on;
% r_tem2=mod.mod{i_mod}.get_r_from_a([0,0,1],1,mod.mod{i_mod}.var.O(:,:,i_c+1));
% quiver3(mod.mod{i_mod}.var.coord(i_c+1,1),mod.mod{i_mod}.var.coord(i_c+1,2),mod.mod{i_mod}.var.coord(i_c+1,3),r_tem2(1),r_tem2(2),r_tem2(3),...
%         10,'linewidth',2);
