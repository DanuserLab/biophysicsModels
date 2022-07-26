classdef ComMath
    methods(Static)
        [G] = grid3(x,y,z,varargin);
        [id_mesh] = getMeshID(coord,Mesh_coord, Mesh_range,Mesh_d,varargin);
        xyzR = rotAxis(r,theta,phi, varargin);
        [j] = getMeshIDNeighbor(i,Mesh_range,varargin);
        [I,check]=plane_line_intersect(n,V0,P0,P1);
        Par = fit2Dcircle(XY);
        [n] = numDigitInt(N,varargin);
    end
end

