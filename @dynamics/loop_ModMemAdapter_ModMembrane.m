function [mod,dyn,abnormal_length] = loop_ModMemAdapter_ModMembrane(dyn,mod,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('dyn', @(x) isa(x,'dynamics'));
            ip.addRequired('mod', @(x) isa(x,'model'));
            ip.addParameter('s_max', 5, @isnumeric);
            ip.addParameter('i_typ', 2, @isnumeric);
            ip.addParameter('i_int', 1, @isnumeric);
            ip.addParameter('ignore_abnormal_length', false, @islogical);
            ip.addParameter('add_spring', false, @islogical);
            ip.addParameter('split', [], @isstruct);
            ip.addParameter('merge', [], @isstruct);
            ip.parse(dyn,mod,varargin{:});
%--------------------------------------------------------------------------
            s_max=ip.Results.s_max;
            i_typ=ip.Results.i_typ;
            i_int=ip.Results.i_int;
            int_name='ModMemAdapter_ModMembrane';
            i_step=dyn.i_step;
            ignore_abnormal_length=ip.Results.ignore_abnormal_length;
%--------------------------------------------------------------------------
%%
            i_mod=[mod.i_mod.ModMembrane,mod.i_mod.ModMemAdapter];  %1: membrane, 2: membrane adapter
%==========================================================================

%%

             [mod.TypForce,loc_relaxed,abnormal_length,mod.mod{i_mod(1)}] = ModMemAdapter_ModMembrane(mod.TypForce,mod.mod{i_mod(1)},mod.mod{i_mod(2)},...
                 'f_const_only',false,'mex',mod.Mex);%,'add_spring',ip.Results.add_spring,'split',ip.Results.split,'merge',ip.Results.merge);

%%
         [mod.TypForce] = var_dt(mod.TypForce,mod.mod{i_mod(1)},'f_const_only',false,'case_f',4);
         [dt_mem,id_tem]=min([mod.TypForce.dt,((floor(mod.t/dyn.pm.dt_const)+1)*dyn.pm.dt_const-mod.t)]);
         dyn.pm.dt=[dyn.pm.dt,dt_mem];
%==========================================================================          
         [dyn] = get_dt(dyn);
         mod.t=mod.t+dyn.pm.dt_var;
%==========================================================================       
%==========================================================================
%%
         mod.mod{i_mod(1)}.var.coord(mod.mod{i_mod(1)}.var.id_on_coord,:)=mod.mod{i_mod(1)}.var.coord(mod.mod{i_mod(1)}.var.id_on_coord,:)...
                                    +(mod.TypForce.int_const.ModMembrane(mod.mod{i_mod(1)}.var.id_on_coord,:)...
                                     +mod.TypForce.int_comp.ModMemAdapter_ModMembrane(mod.mod{i_mod(1)}.var.id_on_coord,:)...
                                     )*mod.mod{i_mod(1)}.pm.mu*dyn.pm.dt_var;
%--------------------------------------------------------------------------                                 
%==========================================================================                                 
         if (abnormal_length==true) && (ignore_abnormal_length==false)
            mod.changed(i_mod(1))=true;
         end
end