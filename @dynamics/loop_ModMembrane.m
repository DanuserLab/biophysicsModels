function [mod,abnormal_length,dyn,V] = loop_ModMembrane(dyn,mod,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('dyn', @(x) isa(x,'dynamics'));
            ip.addRequired('mod', @(x) isa(x,'model'));
            ip.addParameter('f_const_only', false, @islogical);
            ip.addParameter('local', 0, @isnumeric); %0: global; 1: regular local; 2: extension; 3: shrinkage
            ip.addParameter('s_max', 0.01, @isnumeric);
            ip.addParameter('i_typ', 2, @isnumeric);
            ip.addParameter('i_int', 1, @isnumeric);
            ip.addParameter('edg_exo', [], @isnumeric);
            ip.addParameter('ignore_abnormal_length', false, @islogical);
            ip.parse(dyn,mod,varargin{:});
%--------------------------------------------------------------------------
            f_const_only=ip.Results.f_const_only;
            local=ip.Results.local;
            s_max=ip.Results.s_max;
            i_typ=ip.Results.i_typ;
            i_int=ip.Results.i_int;
            ignore_abnormal_length=ip.Results.ignore_abnormal_length;
            edg_exo=ip.Results.edg_exo;
%--------------------------------------------------------------------------
            i_mod=mod.i_mod.ModMembrane;  %membrane
%--------------------------------------------------------------------------
%%
            V=[];
%             for i_step=1:dyn.pm.n_step
                    if dyn.update==true
                       [mod.TypForce,abnormal_length,mod.mod{i_mod},Vtot] = ModMembrane(mod.TypForce,mod,...
                                                          'mex',mod.Mex);  
                       [mod.TypForce,V_tot] = ModMembrane_ModSubstrate(mod.TypForce,mod);
                       dyn.update=false;
                    else
                       [mod.TypForce,abnormal_length,mod.mod{i_mod},Vtot] = ModMembrane(mod.TypForce,mod,...
                                                         'mex',mod.Mex);  
                    end
                                      
                    [mod.TypForce] = var_dt(mod.TypForce,mod.mod{i_mod},'f_const_only',f_const_only,'case_f',1);
                    [dt_tem,id_tem]=min([mod.TypForce.dt,((floor(mod.t/mod.mod{i_mod}.pm.dt)+1)*mod.mod{i_mod}.pm.dt-mod.t)]);
                    if id_tem==2
                        dt_tem=dt_tem+1e-15;
                        [mod.TypForce,abnormal_length,mod.mod{i_mod},Vtot] = ModMembrane(mod.TypForce,mod,...
                                                          'mex',mod.Mex);  
                        [mod.TypForce,V_tot] = ModMembrane_ModSubstrate(mod.TypForce,mod);
                    end
%                     V=[V;V_tot];
   
                    mod.mod{i_mod}.var.coord(mod.mod{i_mod}.var.id_on_coord,:)=mod.mod{i_mod}.var.coord(mod.mod{i_mod}.var.id_on_coord,:)...
                        +(mod.TypForce.int_const.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:)...
                        +mod.TypForce.int_comp.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:)...
                        +mod.TypForce.int_comp.ModMembrane_ModSubstrate{1}(mod.mod{i_mod}.var.id_on_coord,:)...
                        +mod.TypForce.int_rand.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:))*mod.mod{i_mod}.pm.mu*dt_tem;
                    mod.t=mod.t+dt_tem;
dyn.pm.dt=dt_tem;           
m=mod.mod{i_mod};
f=mod.TypForce;
d = (m.var.coord(m.var.edge_all(m.var.id_on_edg,2),:) - m.var.coord(m.var.edge_all(m.var.id_on_edg,1),:));
r = sqrt(sum(d.^2,2));
i_shift=f.int_stored.ModMembrane.rn(1)/m.pm.dr-1;
i = floor(r/m.pm.dr+0.5)-i_shift;
if min(i)<=0
    disp(min(i));
    disp('wall failed');
end
                %fprintf('step: %d --> %d, %f, %d\n',i_step,dyn.pm.n_step,dt_tem,mod.mod{i_mod}.var.n_coord);
%             end
            %fprintf('step: %d ; %d, \n',i_step,mod.mod{i_mod}.var.n_coord);
%             dyn.i_step=dyn.i_step+i_step;    
            dyn.i_step=dyn.i_step+1;
        end