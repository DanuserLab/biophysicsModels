function [M,ftot] = follow(M,nameV,nameIdNeighbor,idModFollower,idModFollowee,ftot,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('M', @(x) isa(x,'model'));
            ip.addRequired('nameV', @(x) ischar(x));
            ip.addRequired('nameIdNeighbor', @(x) ischar(x));
            ip.addRequired('idModFollower', @(x) isnumeric(x));
            ip.addRequired('idModFollowee', @(x) isnumeric(x));
            ip.addRequired('ftot', @(x) iscell(x));
            ip.addParameter('update', false, @islogical);
            ip.addParameter('sMax', 0.01, @isnumeric);
            ip.parse(M,nameV,nameIdNeighbor,idModFollower,idModFollowee,ftot,varargin{:});
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
idCoordFollower=(1:M.mod{idModFollower}.var.n_coord)';
Ner=numel(idCoordFollower);
update=ip.Results.update;
if update==false
    idCoordFollowee=zeros(M.mod{idModFollower}.var.n_coord,1);
else
    idCoordFollowee=M.mod{idModFollower}.var.follow.idCoordFollowee;
end

followExtra=false;
if isa(M.mod{idModFollowee},'ModMembrane')
    idExtra=M.mod{idModFollowee}.var.j_T;
    followExtra=true;
end
% followExtra=false; %XXX

if update==true
    M.mod{idModFollower}.var.coord=M.mod{idModFollowee}.var.coord(idCoordFollowee,:);
    for i=1:M.mod{idModFollower}.var.n_coord
    ftot{idModFollowee}(idCoordFollowee(i),:)=ftot{idModFollowee}(idCoordFollowee(i),:)+ftot{idModFollower}(i,:);
    if followExtra==true
        idExtraAtVer=idExtra(idCoordFollowee(i),:);
        idExtraAtVer=idExtraAtVer(~isnan(idExtraAtVer));
        ftot{idModFollowee}(idExtraAtVer,:)=ftot{idModFollowee}(idExtraAtVer,:)+ftot{idModFollower}(i,:);
    end
    end
else
%see if moving the follower to the neighbors of the followee would lower
%the given potential
for i=1:Ner
%     disp(i);
    d=sum((M.mod{idModFollower}.var.coord(idCoordFollower(i),:)-M.mod{idModFollowee}.var.coord).^2,2);
    [~,M.mod{idModFollower}.var.follow.idCoordFollowee(i)]=min(d);
    M.mod{idModFollower}.var.coord(i,:)=M.mod{idModFollowee}.var.coord(M.mod{idModFollower}.var.follow.idCoordFollowee(i),:);
    ftot{idModFollowee}(M.mod{idModFollower}.var.follow.idCoordFollowee(i),:)=ftot{idModFollowee}(M.mod{idModFollower}.var.follow.idCoordFollowee(i),:)...
                                                                             +ftot{idModFollower}(i,:);
    if followExtra==true
        idExtraAtVer=idExtra(M.mod{idModFollower}.var.follow.idCoordFollowee(i),:);
        idExtraAtVer=idExtraAtVer(~isnan(idExtraAtVer));
        ftot{idModFollowee}(idExtraAtVer,:)=ftot{idModFollowee}(idExtraAtVer,:)+ftot{idModFollower}(i,:);
    end
end
end
[~,V_tot]=M.TypForce.(nameV)(M);
for i=1:Ner
    IdNeighbor=M.mod{idModFollowee}.var.(nameIdNeighbor)(M.mod{idModFollower}.var.follow.idCoordFollowee(i),:);
    nNeighbor=numel(IdNeighbor);
    V_alt=zeros(nNeighbor,1);
    coordTry=zeros(nNeighbor,3);
    coordSave=M.mod{idModFollower}.var.coord(idCoordFollower(i),:);
    for j=1:nNeighbor
        if ~isnan(IdNeighbor(j))
            coordTry(j,:)=M.mod{idModFollowee}.var.coord(IdNeighbor(j),:);
            M.mod{idModFollower}.var.coord(idCoordFollower(i),:)=M.mod{idModFollowee}.var.coord(IdNeighbor(j),:);
            [~,V_alt(j)]=M.TypForce.(nameV)(M,'Vonly',true);            
            M.mod{idModFollower}.var.coord(idCoordFollower(i),:)=coordSave;
        else
            V_alt(j)=inf;
        end
    end
    [Vmin,idMin]=min(V_alt);
    if (Vmin<V_tot)&&(Vmin~=0)
        M.mod{idModFollower}.var.coord(idCoordFollower(i),:)=coordTry(idMin,:);
        M.mod{idModFollower}.var.follow.idCoordFollowee(i)=IdNeighbor(idMin);
    end
end
%--------------------------------------------------------------------------                             
end