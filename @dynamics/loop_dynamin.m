function [mod,dyn,abnormal_length] = loop_dynamin(dyn,mod,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('dyn', @(x) isa(x,'dynamics'));
            ip.addRequired('mod', @(x) isa(x,'model'));
            ip.addParameter('s_max', 5, @isnumeric);
            ip.addParameter('i_typ', 2, @isnumeric);
            ip.addParameter('i_int', 1, @isnumeric);
            ip.addParameter('ignore_abnormal_length', false, @islogical);
            ip.addParameter('add_spring', false, @islogical);
            ip.addParameter('move_adp', true, @islogical);
            ip.addParameter('split', [], @isstruct);
            ip.addParameter('merge', [], @isstruct);
            ip.parse(dyn,mod,varargin{:});
%--------------------------------------------------------------------------
            s_max=ip.Results.s_max;
            i_typ=ip.Results.i_typ;
            i_int=ip.Results.i_int;
            int_name='ModClathrin_ModMemAdapter_ModMembrane';
            i_step=dyn.i_step;
            ignore_abnormal_length=ip.Results.ignore_abnormal_length;
%--------------------------------------------------------------------------
%%
            i_mod=[mod.i_mod.ModClathrin,mod.i_mod.ModMembrane,mod.i_mod.ModMemAdapter];  %1: clathrin, 2: membrane 3:MemAdapter
%==========================================================================
%%
%          if i_step==1
%              [mod.TypForce,loc_relaxed,abnormal_length] = ModMembrane(mod.TypForce,mod.mod{i_mod(2)},...
%                  'f_const_only',false,'mex',mod.Mex);
%              [mod.TypForce] = ModClathrin_ModMemAdapter_ModMembrane(mod.TypForce,mod,i_mod,'ModClathrin_ModMemAdapter_ModMembrane','update',true);
%              [mod.TypForce] = ModClathrin(mod.TypForce,mod,i_mod(1),'ModClathrin');
%          else
%              [mod.TypForce,loc_relaxed,abnormal_length] = ModMembrane(mod.TypForce,mod.mod{i_mod(2)},...
%                  'f_const_only',true,'mex',mod.Mex);
%          end     
if i_step==1
    update_clathrin_membrane=true;
else
    update_clathrin_membrane=false;
end
%%
update_clathrin_membrane=true;

%              [mod.TypForce,loc_relaxed,abnormal_length,mod.mod{i_mod(2)}] = ModMembrane(mod.TypForce,mod.mod{i_mod(2)},...
%                  'f_const_only',false,'mex',mod.Mex);%,'add_spring',ip.Results.add_spring,'split',ip.Results.split,'merge',ip.Results.merge);
             [mod.TypForce,loc_relaxed,abnormal_length,mod.mod{i_mod(2)}] = ModMemAdapter_ModMembrane(mod.TypForce,mod.mod{i_mod(2)},mod.mod{i_mod(3)},...
                 'f_const_only',false,'mex',mod.Mex);

             [mod.TypForce,mod.TypForce.int_V.ModClathrin_ModMemAdapter_ModMembrane]...
                 = ModClathrin_ModMemAdapter_ModMembrane(mod.TypForce,mod,'ModClathrin_ModMemAdapter_ModMembrane','update',update_clathrin_membrane,'mex',mod.Mex);
             [mod.TypForce,mod.TypForce.int_V.ModMembrane_ModSubstrate] = ModMembrane_ModSubstrate(mod.TypForce,mod,'ModMembrane_ModSubstrate');
             [mod.TypForce,mod.TypForce.int_V.ModClathrin] = ModClathrin(mod.TypForce,mod,'ModClathrin','mex',mod.Mex);

%%
% id_feet=[49 50 51];
% for i_test=1:mod.mod{i_mod(1)}.var.n_coord
%     hold on;
% foot=mod.mod{i_mod(1)}.get_r_from_a(mod.mod{i_mod(1)}.var.coord_org(id_feet,1:3),3,mod.mod{i_mod(1)}.var.O(:,:,i_test))...
%      +mod.mod{i_mod(1)}.var.coord(i_test,:);
% f_test=mod.TypForce.int_comp.ModClathrin_ModMemAdapter_ModMembrane{1}(i_test,:);
% quiver3(foot(:,1),foot(:,2),foot(:,3),f_test(:,1),f_test(:,2),f_test(:,3),'linewidth',2); 
% end
%==========================================================================
%%
         [mod.TypForce] = var_dt(mod.TypForce,mod.mod{i_mod(2)},'f_const_only',false,'case_f',6);
         [dt_mem,id_tem]=min([mod.TypForce.dt,((floor(mod.t/dyn.pm.dt_const)+1)*dyn.pm.dt_const-mod.t)]);
         dyn.pm.dt=[dyn.pm.dt,dt_mem];
         
         
%          if id_tem==2
%              [mod.TypForce,loc_relaxed,abnormal_length] = ModMembrane(mod.TypForce,mod.mod{i_mod(2)},...
%                  'f_comp_only',true,'local',local,'mex',mod.Mex);
%              [mod.TypForce] = ModClathrin_ModMemAdapter_ModMembrane(mod.TypForce,mod,'ModClathrin_ModMemAdapter_ModMembrane','update',mod.changed(i_mod(1)));
%              [mod.TypForce] = ModClathrin(mod.TypForce,mod,i_mod(1),'ModClathrin');
%          end
%==========================================================================          
%%
dr=zeros(mod.mod{i_mod(1)}.var.n_coord,6);
da=zeros(mod.mod{i_mod(1)}.var.n_coord,6);
         for ic=1:mod.mod{i_mod(1)}.var.n_coord
             A=mod.mod{i_mod(1)}.var.O(:,:,ic);
             E=mod.mod{i_mod(1)}.var.E(:,:,ic);
             %-------------------------------------------------------------
             mu_r=E*A*mod.mod{i_mod(1)}.pm.mu_r*A'*E';
             mu_r_s=sqrtm(mu_r);
             mu_t=E*A*mod.mod{i_mod(1)}.pm.mu_t*A'*E';
             mu_t_s=sqrtm(mu_t);
             %-------------------------------------------------------------
             Tau=mod.TypForce.int_comp.ModClathrin(ic,4:6)+mod.TypForce.int_comp.ModClathrin_ModMemAdapter_ModMembrane{1}(ic,4:6);
             %a_tem=mod.mod{i_mod(1)}.var.a(ic,:); a_tem=a_tem';
             F=mod.TypForce.int_comp.ModClathrin(ic,1:3)+mod.TypForce.int_comp.ModClathrin_ModMemAdapter_ModMembrane{1}(ic,1:3);
             %r_tem=mod.mod{i_mod(1)}.var.coord(ic,:); r_tem=r_tem';
             da_tem=ones(3,1);dr_tem=ones(3,1);
             Qr=(mu_r*Tau');
             Qt=(mu_t*F');
             dt_tem=dyn.pm.dt_const;
             while (norm(da_tem)>s_max) || (norm(dr_tem)>s_max)
             da_tem=Qr*dt_tem+mu_r_s*randn(3,1)*sqrt(2*mod.pm.kBT*dt_tem);
             dr_tem=(Qt*dt_tem+mu_t_s*randn(3,1)*sqrt(2*mod.pm.kBT*dt_tem))*dyn.pm.cm_to_nm;
%              da_tem=Qr*dt_tem;
%              dr_tem=Qt*dt_tem;
             dt_tem=dt_tem*0.5;
             end
             dyn.pm.dt=[dyn.pm.dt,dt_tem];
             dr(ic,1:3)=Qt';    dr(ic,4:6)=(mu_t_s*randn(3,1))'*0;            
             da(ic,1:3)=Qr';    da(ic,4:6)=(mu_r_s*randn(3,1))'*0;
         end
%========================================================================== 
%%
         [dyn] = get_dt(dyn);
         mod.t=mod.t+dyn.pm.dt_var;
%==========================================================================       
%%
         for ic=1:mod.mod{i_mod(1)}.var.n_coord
             mod.mod{i_mod(1)}.var.a(ic,:)=mod.mod{i_mod(1)}.var.a(ic,:)+da(ic,1:3)*dyn.pm.dt_var+da(ic,4:6)*sqrt(2*1*dyn.pm.dt_var);
             mod.mod{i_mod(1)}.var.ang_a.Phi(ic)=sqrt(sum((mod.mod{i_mod(1)}.var.a(ic,:)).^2,2));
             if mod.mod{i_mod(1)}.var.ang_a.Phi(ic)>pi
                 mod.mod{i_mod(1)}.var.ang_a.Phi(ic)=2*pi-mod.mod{i_mod(1)}.var.ang_a.Phi(ic);
                 mod.mod{i_mod(1)}.var.a(ic,:)=-mod.mod{i_mod(1)}.var.a(ic,:)/norm(mod.mod{i_mod(1)}.var.a(ic,:))*mod.mod{i_mod(1)}.var.ang_a.Phi(ic);
             end
             if (mod.mod{i_mod(1)}.var.ang_a.Phi(ic)>pi) || (mod.mod{i_mod(1)}.var.ang_a.Phi(ic)<0)
                 mod.mod{i_mod(1)}.var.ang_a.Phi(ic)
                 error('wrong a')
             end         
             mod.mod{i_mod(1)}.var.coord(ic,:)=mod.mod{i_mod(1)}.var.coord(ic,:)+dr(ic,1:3)*dyn.pm.dt_var+dr(ic,4:6)*sqrt(2*1*dyn.pm.dt_var);
         end   
         [mod.mod{i_mod(1)}.var] = setVar(mod.mod{i_mod(1)},true);
%==========================================================================
%%
         mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(2)}.var.id_on_coord,:)=mod.mod{i_mod(2)}.var.coord(mod.mod{i_mod(2)}.var.id_on_coord,:)...
                                    +(mod.TypForce.int_const.ModMembrane(mod.mod{i_mod(2)}.var.id_on_coord,:)...
                                     +mod.TypForce.int_comp.ModMemAdapter_ModMembrane(mod.mod{i_mod(2)}.var.id_on_coord,:)...
                                     +mod.TypForce.int_comp.ModClathrin_ModMemAdapter_ModMembrane{2}(mod.mod{i_mod(2)}.var.id_on_coord,:) ...
                                     +mod.TypForce.int_comp.ModMembrane_ModSubstrate{1}(mod.mod{i_mod(2)}.var.id_on_coord,:) ...
                                     )*mod.mod{i_mod(2)}.pm.mu*dyn.pm.dt_var;
%--------------------------------------------------------------------------         
         if ip.Results.move_adp==true
         id_m_all=(1:mod.mod{i_mod(2)}.var.n_coord)';
         id_m_for_adp_change=id_m_all(mod.TypForce.int_comp.(int_name){4}>0);
         n_adp_change=numel(id_m_for_adp_change);
         if n_adp_change>0   
         id_ModMembrane_new=mod.mod{i_mod(3)}.var.id_ModMembrane;
         occupied_new=mod.mod{i_mod(3)}.var.occupied;
         for i=1:n_adp_change
             if ismember(mod.TypForce.int_comp.(int_name){4}(id_m_for_adp_change(i)),mod.mod{i_mod(2)}.var.id_on_coord)
                   ardy_adp=ismember(mod.TypForce.int_comp.(int_name){4}(id_m_for_adp_change(i)),id_ModMembrane_new);
                   if  (ardy_adp==false)
                       id_adp_change=mod.mod{i_mod(3)}.var.id_ModMembrane==id_m_for_adp_change(i);
                       %id_adp_change=id_adp_change.*(~mod.mod{i_mod(3)}.var.occupied);
                       id_adp_change=id_adp_change==1;
                       id_ModMembrane_new(id_adp_change)=mod.TypForce.int_comp.(int_name){4}(id_m_for_adp_change(i));
%                        id_adp_new=id_ModMembrane_new==mod.TypForce.int_comp.(int_name){4}(id_m_for_adp_change(i));
%                        id_adp_new=id_adp_new==1;
                       %occupied_new(id_adp_change)=false;
                       %occupied_new(id_adp_new)=true;
                   else
                       %id_adp_new=mod.mod{i_mod(3)}.var.id_ModMembrane==mod.TypForce.int_comp.(int_name){4}(id_m_for_adp_change(i));
                       id_adp_new=id_ModMembrane_new==mod.TypForce.int_comp.(int_name){4}(id_m_for_adp_change(i));
                       id_adp_new=id_adp_new==1;
                       if mod.mod{i_mod(3)}.var.occupied(id_adp_new)==false
                           id_adp_change=mod.mod{i_mod(3)}.var.id_ModMembrane==id_m_for_adp_change(i);
                           %id_adp_change=id_adp_change.*(~mod.mod{i_mod(3)}.var.occupied);
                           id_adp_change=id_adp_change==1;
                           id_ModMembrane_new(id_adp_change)=mod.TypForce.int_comp.(int_name){4}(id_m_for_adp_change(i));
                           id_ModMembrane_new(id_adp_new)=id_m_for_adp_change(i);
                           %occupied_new(id_adp_change)=false;
                           %occupied_new(id_adp_new)=true;
                       end
                   end
             end
         end
         mod.mod{i_mod(3)}.var.id_ModMembrane=id_ModMembrane_new; 
         mod.mod{i_mod(3)}.var.occupied=occupied_new;
         end
         end
%==========================================================================                                 
         if (abnormal_length==true) && (ignore_abnormal_length==false)
            mod.changed(i_mod(2))=true;
         end
end