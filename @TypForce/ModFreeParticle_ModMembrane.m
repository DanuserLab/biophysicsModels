function [f,V_tot] = ModFreeParticle_ModMembrane(f,mod,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addParameter('idClathrinSub', [], @isnumeric);
ip.parse(f,mod,varargin{:}); 
%----------------------------------------------------------------------------------------
            i_mod=[mod.i_mod.ModFreeParticle,mod.i_mod.ModMembrane];  %1: AP2 as free particle; 2: ModMembrane
%----------------------------------------------------------------------------------------
mod.mod{mod.i_mod.ModFreeParticle}.var.idMesh = ...
ComMath.getMeshID(mod.mod{mod.i_mod.ModFreeParticle}.var.coord,mod.Mesh.coord, mod.Mesh.range,mod.Mesh.d);
mod.mod{mod.i_mod.ModMembrane}.var.idMesh = ...
ComMath.getMeshID(mod.mod{mod.i_mod.ModMembrane}.var.coord,mod.Mesh.coord, mod.Mesh.range,mod.Mesh.d);
%----------------------------------------------------------------------------------------
idMeshFreeParticle=mod.mod{mod.i_mod.ModFreeParticle}.var.idMesh';

[idMeshMembrane] = ComMath.getMeshIDNeighbor(mod.mod{mod.i_mod.ModMembrane}.var.idMesh',mod.Mesh.range);

idMembrane=[];
idFreeParticle=[];
for i=1:27 %is nMeshFreeParticle=size(idMeshFreeParticle,1); 27 neighbors at one mesh point
IdM=(idMeshFreeParticle-idMeshMembrane(i,:)');

[idMembraneTem,idFreeParticleTem]=find(IdM==0);
idMembrane=[idMembrane;idMembraneTem];
idFreeParticle=[idFreeParticle;idFreeParticleTem];
end

f.int_comp.ModFreeParticle_ModMembrane=cell(2,1);
f.int_comp.ModFreeParticle_ModMembrane{1}=zeros(mod.mod{i_mod(1)}.var.n_coord,3);
f.int_comp.ModFreeParticle_ModMembrane{2}=zeros(mod.mod{i_mod(2)}.var.n_coord,3);

nPair=numel(idMembrane);
Vpair=zeros(nPair,1);
Fpair=zeros(nPair,3);

for iPair=1:nPair
   Vpair(iPair)=0.5*f.pm.k_ModFreeParticle_ModMembrane*...
          sum((mod.mod{mod.i_mod.ModMembrane}.var.coord(idMembrane(iPair),:)-mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle(iPair),:)).^2,2);
   Fpair(iPair,:) = -f.pm.k_ModFreeParticle_ModMembrane*(mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle(iPair),:)-...
                                                       mod.mod{mod.i_mod.ModMembrane}.var.coord(idMembrane(iPair),:));%force on ModFreeParticle 
end

[Cunq,~,iCunq]=unique(idFreeParticle,'row');

Nunq=size(Cunq,1);

Vunq=inf(Nunq,1);

idMin=zeros(Nunq,1);

for iPair=1:nPair
    if Vunq(iCunq(iPair)) > Vpair(iPair)
        Vunq(iCunq(iPair))=Vpair(iPair);
        idMin(iCunq(iPair))=iPair;
    end
end
for iUnq=1:Nunq
    iPair=idMin(iUnq);
    f.int_comp.ModFreeParticle_ModMembrane{1}(idFreeParticle(iPair),:)=Fpair(iPair,:);
    f.int_comp.ModFreeParticle_ModMembrane{2}(idMembrane(iPair),:)=-Fpair(iPair,:);
end
%taking 1-ring average
f_tem=zeros(size(f.int_comp.ModFreeParticle_ModMembrane{2}));
for i=1:mod.mod{mod.i_mod.ModMembrane}.var.n_coord
    n_node=mod.mod{mod.i_mod.ModMembrane}.var.n_node(i);
    f_tem(i,:)=mean(f.int_comp.ModFreeParticle_ModMembrane{2}(mod.mod{mod.i_mod.ModMembrane}.var.j_T(i,1:n_node),:));
end
f.int_comp.ModFreeParticle_ModMembrane{2}=f_tem;
V_tot=sum(Vunq);
f.int_V.ModFreeParticle_ModMembrane=V_tot;
f.int_tot.ModFreeParticle_ModMembrane=f.int_comp.ModFreeParticle_ModMembrane;
%--------------------------------------------------------------------------
%%
% figure();
% quiver3(mod.mod{mod.i_mod.ModMembrane}.var.coord(:,1),...
%         mod.mod{mod.i_mod.ModMembrane}.var.coord(:,2),...
%         mod.mod{mod.i_mod.ModMembrane}.var.coord(:,3),...
%         f.int_tot.ModFreeParticle_ModMembrane{2}(:,1),...
%         f.int_tot.ModFreeParticle_ModMembrane{2}(:,2),...
%         f.int_tot.ModFreeParticle_ModMembrane{2}(:,3),'linewidth',2);
% fig=figure('units','normalized','outerposition',[0 0 1 1]);
%             plot(mod.mod{mod.i_mod.ModMembrane},'linestyle','-','f',fig,'FaceAlpha',1);
%             scatter3(mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle,1),...
%                      mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle,2),...
%                      mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle,3),10,'filled');
%             hold on;
%             scatter3(mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle(idMin),1),...
%                      mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle(idMin),2),...
%                      mod.mod{mod.i_mod.ModFreeParticle}.var.coord(idFreeParticle(idMin),3),20,'filled');
%                  x_lim=[-15 15];
%            xlim(x_lim);ylim(x_lim);zlim(x_lim);
%--------------------------------------------------------------------------
end


