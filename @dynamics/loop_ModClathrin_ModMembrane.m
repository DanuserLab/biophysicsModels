function [mod] = loop_ModClathrin_ModMembrane(dyn,mod,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('dyn', @(x) isa(x,'dynamics'));
            ip.addRequired('mod', @(x) isa(x,'model'));
            ip.addParameter('s_max', 0.01, @isnumeric);
            ip.addParameter('i_typ', 2, @isnumeric);
            ip.addParameter('i_int', 1, @isnumeric);
            ip.addParameter('i_mod', [1 2], @isnumeric);
            ip.parse(dyn,mod,varargin{:});
%--------------------------------------------------------------------------
            s_max=ip.Results.s_max;
            i_typ=ip.Results.i_typ;
            i_int=ip.Results.i_int;
            i_mod=ip.Results.i_mod;
            int_name='ModClathrin_ModMembrane';
            i_step=dyn.i_step;
%--------------------------------------------------------------------------
%x_save=zeros(dyn.pm.n_step,4);

%          mod.mod{i_mod(1)}.var.E=mod.mod{i_mod(1)}.getE(mod.mod{i_mod(1)}.var.a,mod.mod{i_mod(1)}.var.ang_a.Phi);
%          mod.mod{i_mod(1)}.var.O=mod.mod{i_mod(1)}.Omega(mod.mod{i_mod(1)}.var.a,mod.mod{i_mod(1)}.var.ang_a.Phi);
         [mod.TypForce] = ModClathrin_ModMembrane(mod.TypForce,mod,i_mod,'ModClathrin_ModMembrane','update',mod.changed(i_mod(1)));
         [mod.TypForce] = ModClathrin(mod.TypForce,mod,i_mod(1),'ModClathrin');
         for ic=1:mod.mod{i_mod(1)}.var.n_coord
             A=mod.mod{i_mod(1)}.var.O(:,:,ic);
             E=mod.mod{i_mod(1)}.var.E(:,:,ic);
             %-------------------------------------------------------------
             mu_r=E*A*mod.mod{i_mod(1)}.pm.mu_r*A'*E';
             mu_r_s=sqrtm(mu_r);
             mu_t=E*A*mod.mod{i_mod(1)}.pm.mu_t*A'*E';
             mu_t_s=sqrtm(mu_t);
             %-------------------------------------------------------------
             Tau=mod.TypForce.int_comp.ModClathrin(ic,4:6)+mod.TypForce.int_comp.ModClathrin_ModMembrane{1}(ic,4:6);
             a_tem=mod.mod{i_mod(1)}.var.a(ic,:);
             a_tem=a_tem';
             F=mod.TypForce.int_comp.ModClathrin(ic,1:3)+mod.TypForce.int_comp.ModClathrin_ModMembrane{1}(ic,1:3);
             r_tem=mod.mod{i_mod(1)}.var.coord(ic,:);
             r_tem=r_tem';
             dt_tem=dyn.pm.dt;
             da=ones(3,1);dr=ones(3,1);
             Qr=(mu_r*Tau');
             Qt=(mu_t*F');
             while (norm(da)>s_max) || (norm(dr)>s_max)
             %da=Qr*dt_tem+mu_r_s*randn(3,1)*sqrt(2*mod.pm.kBT*dt_tem);
             %dr=(Qt*dt_tem+mu_t_s*randn(3,1)*sqrt(2*mod.pm.kBT*dt_tem))*dyn.pm.cm_to_nm;
             da=Qr*dt_tem;
             dr=Qt*dt_tem*dyn.pm.cm_to_nm;
             dt_tem=dt_tem*0.5;
             end
             dyn.dt=[dyn.dt;dt_tem];
             r_tem=r_tem+dr;            
             a_tem=a_tem+da;
             mod.mod{i_mod(1)}.var.a(ic,:)=a_tem';
             mod.mod{i_mod(1)}.var.ang_a.Phi(ic)=sqrt(sum((mod.mod{i_mod(1)}.var.a(ic,:)).^2,2));
             if mod.mod{i_mod(1)}.var.ang_a.Phi(ic)>pi
                 mod.mod{i_mod(1)}.var.ang_a.Phi(ic)=2*pi-mod.mod{i_mod(1)}.var.ang_a.Phi(ic);
                 mod.mod{i_mod(1)}.var.a(ic,:)=-mod.mod{i_mod(1)}.var.a(ic,:)/norm(mod.mod{i_mod(1)}.var.a(ic,:))*mod.mod{i_mod(1)}.var.ang_a.Phi(ic);
             end
             if (mod.mod{i_mod(1)}.var.ang_a.Phi(ic)>pi) || (mod.mod{i_mod(1)}.var.ang_a.Phi(ic)<0)
                 mod.mod{i_mod(1)}.var.ang_a.Phi(ic)
                 error('wrong a')
             end         
             mod.mod{i_mod(1)}.var.coord(ic,:)=r_tem';

            %mod.mod{i_mod}.var.O(:,:,ic)=mod.mod{i_mod}.Omega(mod.mod{i_mod}.var.a(ic,:),mod.mod{i_mod}.var.ang_a.Phi(ic));
            %mod.mod{i_mod}.var.E(:,:,ic)=mod.mod{i_mod}.getE(mod.mod{i_mod}.var.a(ic,:),mod.mod{i_mod}.var.ang_a.Phi(ic));
         end
         %-----------------------------------------------------------------
         mod.mod{i_mod(2)}.var.coord=mod.mod{i_mod(2)}.var.coord+1e10*mod.TypForce.int_comp.ModClathrin_ModMembrane{2}*dt_save;
         %-----------------------------------------------------------------

     [mod.mod{i_mod(1)}.var] = setVar(mod.mod{i_mod(1)},true);
%      x_save(i_step,1:3)=mod.mod{i_mod(1)}.var.a(1,:);
%      x_save(i_step,4)=sum(sum(V));
%      figure;
%      subplot(1,2,1);
%      plot(x_save(:,1));hold on;plot(x_save(:,2));hold on;plot(x_save(:,3));
%      subplot(1,2,2);plot(x_save(:,4));
end