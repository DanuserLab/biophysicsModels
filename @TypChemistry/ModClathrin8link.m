function [mod] = ModClathrin8link(ch,mod,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('ch', @(x) isa(x,'TypChemistry'));
ip.addRequired('mod', @(x) isa(x,'model'));
ip.addParameter('mex', [], @isobject);
ip.addParameter('r_min', [], @isnumeric);
ip.parse(ch,mod,varargin{:});
%--------------------------------------------------------------------------
            i_mod=mod.i_mod.ModClathrin;  %1: clathrin
%--------------------------------------------------------------------------
r_min=ip.Results.r_min;
if isempty(r_min)
%  r_min=norm(mod.mod{1}.var.coord_org(1,1:3)-mod.mod{1}.var.coord_org(2,1:3))*10;
   r_min=mod.mod{i_mod}.pm.r_min_for_link;
end
%----------------------------------------------------------------------------------------
%mex_avail=ip.Results.mex_avail;
mex=ip.Results.mex;
if isempty(mex)
    mex_avail=false;
else
    mex_avail=true;
end
%----------------------------------------------------------------------------------------
[waist] = getWaist(mod.mod{i_mod});
% [CtrlPt] = getCtrlPt(mod.mod{i_mod});
% u=zeros(mod.mod{i_mod}.var.n_coord,3);
% for i_c=1:mod.mod{i_mod}.var.n_coord
% u(i_c,:)=mod.mod{i_mod}.get_r_from_a([0,0,1],1,mod.mod{i_mod}.var.O(:,:,i_c));
% end
%%
id_all=(1:3);
for i_c1=1:mod.mod{i_mod}.var.n_coord-1
%     disp(i_c1);
    if ismember(0,mod.mod{i_mod}.var.connect(:,1,i_c1))
    for i_c2=i_c1+1:mod.mod{i_mod}.var.n_coord
%         disp(i_c2);
        %cos_tem=dot(u(i_c1,:),u(i_c2,:));%norm(u(i_c1,:)) and norm(u(i_c2,:)) = 1
        if (~ismember(i_c2,mod.mod{i_mod}.var.connect(:,1,i_c1))) ...
          &&(ismember(0,mod.mod{i_mod}.var.connect(:,1,i_c2)))
%--------------------------------------------------------------------------  
            id1=id_all(mod.mod{i_mod}.var.connect(:,1,i_c1)==0);
            id2=id_all(mod.mod{i_mod}.var.connect(:,1,i_c2)==0);
            d1=sqrt(sum((waist(id1,:,i_c1)-mod.mod{i_mod}.var.coord(i_c2,:)).^2,2));
            d2=sqrt(sum((waist(id2,:,i_c2)-mod.mod{i_mod}.var.coord(i_c1,:)).^2,2));
            id1=id1(d1<r_min);
            id2=id2(d2<r_min);
            if (~isempty(id1)) && (~isempty(id2))
                %d_Ctrl=sqrt(sum((CtrlPt(id1,:,i_c1)-CtrlPt(id2,:,i_c2)).^2,2));
                %if d_Ctrl<r_min
%                 disp('yeah');pause;
                connect_save=mod.mod{i_mod}.var.connect;
                mod.mod{i_mod}.var.connect(id1(1),1,i_c1)=i_c2;
                mod.mod{i_mod}.var.connect(id1(1),2,i_c1)=id2(1);
                mod.mod{i_mod}.var.connect(id2(1),1,i_c2)=i_c1;
                mod.mod{i_mod}.var.connect(id2(1),2,i_c2)=id1(1);
                junction=[i_c1 id1(1);i_c2 id2(1)];
%                 disp('before');
%                 disp(junction);
%                 pause;
                [topology,n_ring] = getTopology(mod.mod{i_mod},junction);
%                 disp('after');
                if (topology==true) && (n_ring>4)
                    fprintf('new connection added: %d,%d; %d,%d\n',i_c1,id1(1),i_c2,id2(1));
                    mod.changed(i_mod)=true;
                else
                    mod.mod{i_mod}.var.connect=connect_save;
                end
                
                %end
            end
%             d=sqrt(sum((CtrlPt(id1,:,i_c1)-CtrlPt(id2,:,i_c2)).^2,2));
%             if d<r_min
%                 mod.mod{i_mod}.var.connect(id1(1),1,i_c1)=i_c2;
%                 mod.mod{i_mod}.var.connect(id1(1),2,i_c1)=id2(1);
%                 mod.mod{i_mod}.var.connect(id2(1),1,i_c2)=i_c1;
%                 mod.mod{i_mod}.var.connect(id2(1),2,i_c2)=id1(1);
%             end
%--------------------------------------------------------------------------
        end
    end
    end
end

