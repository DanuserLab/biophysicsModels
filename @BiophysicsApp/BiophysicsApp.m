classdef BiophysicsApp
    methods(Static)
        [] = TestCMEforce(mod,varargin);
        [] = tetherPull(dirRoot,dir_mod,varargin);
        [] = redBloodCell(dirRoot,dirMod,varargin);
        [] = membraneFussion(dirRoot,dirMod,varargin);
        [] = exampleFilopodia(dirRoot,dirMod,varargin);
        [] = exampleLamellipodia(dirRoot,dirMod,varargin);
        [] = exampleEndocytosis(dirRoot,dirMod,varargin);
    end
end

