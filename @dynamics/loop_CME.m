function [mod,dyn,abnormal_length] = loop_CME(dyn,mod,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('dyn', @(x) isa(x,'dynamics'));
            ip.addRequired('mod', @(x) isa(x,'model'));
            ip.addParameter('add_spring', false, @islogical);
            ip.addParameter('move_adp', true, @islogical);
            ip.addParameter('split', [], @isstruct);
            ip.addParameter('merge', [], @isstruct);
            ip.addParameter('sMax', [], @isnumeric);
            ip.addParameter('sMaxCutoff', 1e-10, @isnumeric);
            ip.addParameter('step_max_mechanics', 1000, @isnumeric);
            ip.parse(dyn,mod,varargin{:});
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
            i_mod=[mod.i_mod.ModClathrin,mod.i_mod.ModMembrane,mod.i_mod.ModFreeParticle];  %1: clathrin, 2: membrane 3:FreeParticle
%--------------------------------------------------------------------------
            sMax=ip.Results.sMax;
            if isempty(sMax)
            sMax=mod.mod{i_mod(2)}.pm.l0*0.1;
            end
            sMaxCutoff=ip.Results.sMaxCutoff;
            step_max_mechanics=ip.Results.step_max_mechanics;
%==========================================================================
%%
cut_off_mechanics=false;
dyn.i_step=0;
VtotPre=inf;
n_AP2=mod.mod{mod.i_mod.ModFreeParticle}.var.n_coord;

while cut_off_mechanics==false
    [mod.TypForce,abnormal_length,mod.mod{i_mod(2)},Vtot] = ModMembrane(mod.TypForce,mod.mod{i_mod(2)},'mex',mod.Mex);
    if abnormal_length==true
%--------------------------------------------------------------------------
        %         fprintf('remesh needed at %d\n',i_par);
        [mod,id_s,id_m] = ModMembrane(mod.TypChemistry,mod,'f_const_only',true,'print_or_not',false);
        rTem2=sum((mod.mod{mod.i_mod.ModMembrane}.var.coord(:,1).^2+mod.mod{mod.i_mod.ModMembrane}.var.coord(:,2).^2),2);
        [~,idTem]=sort(rTem2);
        [mod.mod{mod.i_mod.ModFreeParticle}] = SetVar(mod.mod{mod.i_mod.ModFreeParticle},false,...
            'coordAssigned',mod.mod{mod.i_mod.ModMembrane}.var.coord(idTem(1:n_AP2),:));
    else
%--------------------------------------------------------------------------
        [mod.TypForce,mod.TypForce.int_V.ModClathrin] = ModClathrin(mod.TypForce,mod,'mex',mod.Mex);
        [mod.TypForce,mod.TypForce.int_V.ModClathrin_ModFreeParticle,identifier,idConnectC_FP] = ModClathrin_ModFreeParticle(mod.TypForce,mod);
        Vtot=Vtot+mod.TypForce.int_V.ModClathrin+mod.TypForce.int_V.ModClathrin_ModFreeParticle;
        % [mod.TypForce,V_tot] = ModMembrane_ModSubstrate(mod.TypForce,mod,'ModMembrane_ModSubstrate');
        if Vtot<VtotPre
            VtotPre=Vtot;
            [~,idConnectF_M] = follow(mod.mod{i_mod(3)},mod,'ModMembrane');
            
            idM=idConnectF_M(idConnectC_FP(:,3),2);
            
            fTotonMembrane=zeros(size(mod.mod{i_mod(2)}.var.coord));
            fTotonMembrane(idM,:)=fTotonMembrane(idM,:)+...
                mod.TypForce.int_comp.ModClathrin_ModFreeParticle{identifier(mod.i_mod.ModFreeParticle)}(idConnectC_FP(:,3),:);
            fTotonMembrane=fTotonMembrane+mod.TypForce.int_const.ModMembrane+mod.TypForce.int_comp.ModMembrane;
            
            dt=varDt_ModMembrane(dyn,mod.TypForce,fTotonMembrane,mod.mod{i_mod(2)});
            [coordAttempt,dt_new] = dynamics.translation(mod.mod{i_mod(2)}.var.coord,fTotonMembrane,mod.mod{i_mod(2)}.pm.mu,dt);
            if dt_new<dt
                dt=dt_new;
            end
            [mod.mod{mod.i_mod.ModClathrin},dt_new] = dynamics.rotation(mod.mod{mod.i_mod.ModClathrin},mod.TypForce.int_comp.ModClathrin{1},dt,0,'sMax',sMax);
            if dt_new<dt
                dt=dt_new;
                [mod.mod{i_mod(2)}.var.coord,~] = dynamics.translation(mod.mod{i_mod(2)}.var.coord,fTotonMembrane,mod.mod{i_mod(2)}.pm.mu,dt);
            else
                mod.TypForce,mod.mod{i_mod(2)}.var.coord=coordAttempt;
            end
            rTem2=sum((mod.mod{mod.i_mod.ModMembrane}.var.coord(:,1).^2+mod.mod{mod.i_mod.ModMembrane}.var.coord(:,2).^2),2);
            [~,idTem]=sort(rTem2);
            [mod.mod{mod.i_mod.ModFreeParticle}] = SetVar(mod.mod{mod.i_mod.ModFreeParticle},false,...
                'coordAssigned',mod.mod{mod.i_mod.ModMembrane}.var.coord(idTem(1:n_AP2),:));
            cPre=mod.mod{mod.i_mod.ModClathrin};
            McoordPre=mod.mod{i_mod(2)}.var.coord;
            fpPre=mod.mod{mod.i_mod.ModFreeParticle};
%--------------------------------------------------------------------------
        else
            mod.mod{mod.i_mod.ModClathrin}=cPre;
            mod.mod{i_mod(2)}.var.coord=McoordPre;
            mod.mod{mod.i_mod.ModFreeParticle}=fpPre;
            sMax=sMax*0.1;
%--------------------------------------------------------------------------
        end
    end
    fprintf('%d   %f\n',dyn.i_step,VtotPre);
    if (dyn.i_step > step_max_mechanics) || (sMax<sMaxCutoff)
        cut_off_mechanics=true;
    end
end
%--------------------------------------------------------------------------
end