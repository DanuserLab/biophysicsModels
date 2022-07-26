function [mod,loc_relaxed,abnormal_length] = loop_ModMembrane_bound(dyn,mod,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('dyn', @(x) isa(x,'dynamics'));
            ip.addRequired('mod', @(x) isa(x,'model'));
            ip.addParameter('f_const_only', false, @islogical);
            ip.addParameter('local', 0, @isnumeric); %0: global; 1: regular local; 2: extension; 3: shrinkage
            ip.addParameter('s_max', 0.01, @isnumeric);
            ip.addParameter('i_typ', 2, @isnumeric);
            ip.addParameter('i_int', 1, @isnumeric);
            ip.addParameter('i_mod', 1, @isnumeric);
            ip.addParameter('ignore_abnormal_length', false, @islogical);
            ip.parse(dyn,mod,varargin{:});
%--------------------------------------------------------------------------
            f_const_only=ip.Results.f_const_only;
            local=ip.Results.local;
            s_max=ip.Results.s_max;
            i_typ=ip.Results.i_typ;
            i_int=ip.Results.i_int;
            i_mod=ip.Results.i_mod;
            ignore_abnormal_length=ip.Results.ignore_abnormal_length;
%--------------------------------------------------------------------------
            for i_step=1:dyn.pm.n_step
                if (local==0)
                    if i_step==1
                       [mod.TypForce,loc_relaxed,abnormal_length] = ModMembrane(mod.TypForce,mod.mod{i_mod},...
                                                          'f_const_only',false,'local',local,'mex',mod.Mex);  
                       [mod.TypForce] = ModMembrane_ModSubstrate(mod.TypForce,mod,mod.mat_int{i_typ}{i_int}.i_mod,mod.mat_int{i_typ}{i_int}.name);
                    else
                       [mod.TypForce,loc_relaxed,abnormal_length] = ModMembrane(mod.TypForce,mod.mod{i_mod},...
                                                          'f_const_only',true,'local',local,'mex',mod.Mex);  
                    end
                                      
                    [mod.TypForce] = var_dt(mod.TypForce,mod.mod{i_mod},'f_const_only',f_const_only);
                    [dt_tem,id_tem]=min([mod.TypForce.dt,((floor(mod.t/dyn.pm.dt)+1)*dyn.pm.dt-mod.t)]);
                    if id_tem==2
                        [mod.TypForce,loc_relaxed,abnormal_length] = ModMembrane(mod.TypForce,mod.mod{i_mod},...
                                                          'f_comp_only',true,'local',local,'mex',mod.Mex);  
                        [mod.TypForce] = ModMembrane_ModSubstrate(mod.TypForce,mod,mod.mat_int{i_typ}{i_int}.i_mod,mod.mat_int{i_typ}{i_int}.name);
                    end
               
                    mod.mod{i_mod}.var.coord(mod.mod{i_mod}.var.id_on_coord,:)=mod.mod{i_mod}.var.coord(mod.mod{i_mod}.var.id_on_coord,:)...
                                    +(mod.TypForce.int_const.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:)...
                                     +mod.TypForce.int_comp.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:)...
                                     +mod.TypForce.int_comp.ModMembrane_ModSubstrate{1}(mod.mod{i_mod}.var.id_on_coord,:))*dt_tem;
                else %for local
                    if local==1
                        [mod.TypForce,loc_relaxed,abnormal_length] = ModMembrane(mod.TypForce,mod.mod{i_mod},...
                                                          'f_const_only',f_const_only,'local',local,'mex',mod.Mex);  
                    elseif (local==2) || (local==3)
                        id_tem=sum((mod.mod{i_mod}.var.edge_all(mod.mod{i_mod}.var.id_on_edg,:)==mod.mod{i_mod}.var.id_on_coord(1))+...
                        (mod.mod{i_mod}.var.edge_all(mod.mod{i_mod}.var.id_on_edg,:)==mod.mod{i_mod}.var.id_on_coord(2)),2)==2;
                        edg_exo=mod.mod{i_mod}.var.id_on_edg(id_tem);
                        [mod.TypForce,loc_relaxed,abnormal_length] = ModMembrane(mod.TypForce,mod.mod{i_mod},...
                                                          'f_const_only',f_const_only,'local',local,'mex',mod.Mex,'edg_exo',edg_exo); 
                    end
                    
                    dt_tem=dyn.pm.dt;                   
                    f_mag = max([sqrt(sum(mod.TypForce.int_const.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:).^2,2));...
                        sqrt(sum(mod.TypForce.int_const.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:).^2,2))]);
                    if f_mag*dt_tem > s_max
                        dt_tem = s_max/f_mag;
                    end
                    mod.mod{i_mod}.var.coord(mod.mod{i_mod}.var.id_on_coord,:)=mod.mod{i_mod}.var.coord(mod.mod{i_mod}.var.id_on_coord,:)+mod.TypForce.int_const.ModMembrane(mod.mod{i_mod}.var.id_on_coord,:)*dt_tem;
                end
                if (local==0) && (abnormal_length==true) && (ignore_abnormal_length==false)
                    break;
                elseif (local==1) && (loc_relaxed==true)
                    break;
                elseif (local==2) && ((loc_relaxed==1) || (loc_relaxed==2))
                    break;
                elseif (local==3) && ((loc_relaxed==1) || (loc_relaxed==2))
                    break;
                end               
                %fprintf('step: %d --> %d, %f, %d\n',i_step,dyn.pm.n_step,dt_tem,mod.mod{i_mod}.var.n_coord);
            end
            %fprintf('step: %d ; %d, \n',i_step,mod.mod{i_mod}.var.n_coord);
        end