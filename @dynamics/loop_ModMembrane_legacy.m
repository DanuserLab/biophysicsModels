function [mod,loc_relaxed,abnormal_length,dyn,V] = loop_ModMembrane(dyn,mod,varargin)
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
            for i_step=1:dyn.pm.n_step
                if (local==0)
                    if i_step==1
                       [mod.TypForce,loc_relaxed,abnormal_length,mod.mod{i_mod}] = ModMembrane(mod.TypForce,mod.mod{i_mod},...
                                                          'f_const_only',false,'local',local,'mex',mod.Mex);  
                       [mod.TypForce,V_tot] = ModMembrane_ModSubstrate(mod.TypForce,mod,'ModMembrane_ModSubstrate');
                    else
                       [mod.TypForce,loc_relaxed,abnormal_length] = ModMembrane(mod.TypForce,mod.mod{i_mod},...
                                                          'f_const_only',true,'local',local,'mex',mod.Mex);  
                    end
                                      
                    [mod.TypForce] = var_dt(mod.TypForce,mod.mod{i_mod},'f_const_only',f_const_only,'case_f',1);
                    [dt_tem,id_tem]=min([mod.TypForce.dt,((floor(mod.t/dyn.pm.dt_const)+1)*dyn.pm.dt_const-mod.t)]);
                    if id_tem==2
                        
                        [mod.TypForce,loc_relaxed,abnormal_length,mod.mod{i_mod}] = ModMembrane(mod.TypForce,mod.mod{i_mod},...
                                                          'f_comp_only',true,'local',local,'mex',mod.Mex);  
                        [mod.TypForce,V_tot] = ModMembrane_ModSubstrate(mod.TypForce,mod,'ModMembrane_ModSubstrate');
                    end
                    V=[V;V_tot];
                else %for local
                    if local==1
                        [mod.TypForce,loc_relaxed,abnormal_length] = ModMembrane(mod.TypForce,mod.mod{i_mod},...
                                                          'f_const_only',f_const_only,'local',local,'mex',mod.Mex);  
                    elseif (local==2) || (local==3)
                        [mod.TypForce,loc_relaxed,abnormal_length] = ModMembrane(mod.TypForce,mod.mod{i_mod},...
                                                          'f_const_only',f_const_only,'local',local,'mex',mod.Mex,'edg_exo',edg_exo); 
                    end                   
                end
                if (local==0) && (abnormal_length==true) && (ignore_abnormal_length==false)
                    break;
                elseif (local==1) && (loc_relaxed==true)
                    break;
                elseif (local==2) && ((loc_relaxed==1) || (loc_relaxed==2))
                %elseif (local==2) && ((loc_relaxed==2))
                    break;
                elseif (local==3) && ((loc_relaxed==1) || (loc_relaxed==2))
                %elseif (local==3) && ((loc_relaxed==2))
                    break;
                end    
                if (local==0)
                    mod.mod{i_mod}.var.coord(mod.mod{i_mod}.var.id_on_coord,:)=mod.mod{i_mod}.var.coord(mod.mod{i_mod}.var.id_on_coord,:)...
                        +(mod.TypForce.int_const.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:)...
                        +mod.TypForce.int_comp.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:)...
                        +mod.TypForce.int_comp.ModMembrane_ModSubstrate{1}(mod.mod{i_mod}.var.id_on_coord,:)...
                        +mod.TypForce.int_rand.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:))*mod.mod{i_mod}.pm.mu*dt_tem;
                    mod.t=mod.t+dt_tem;
                else
                    dt_tem=dyn.pm.dt_const;                   
                    f_mag = max([sqrt(sum(mod.TypForce.int_const.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:).^2,2));...
                        sqrt(sum(mod.TypForce.int_const.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:).^2,2))]);
                    if f_mag*dt_tem > s_max
                        dt_tem = s_max/f_mag;
                    end
                    mod.mod{i_mod}.var.coord(mod.mod{i_mod}.var.id_on_coord,:)=mod.mod{i_mod}.var.coord(mod.mod{i_mod}.var.id_on_coord,:)+mod.TypForce.int_const.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:)*dt_tem;
                end
                %fprintf('step: %d --> %d, %f, %d\n',i_step,dyn.pm.n_step,dt_tem,mod.mod{i_mod}.var.n_coord);
            end
            %fprintf('step: %d ; %d, \n',i_step,mod.mod{i_mod}.var.n_coord);
            dyn.i_step=dyn.i_step+i_step;       
        end