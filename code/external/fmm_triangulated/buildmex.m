% mex -output fastmarch1_mex -I./one2all one2all/my_heap.c one2all/unfold.c one2all/marchAlg.c one2all/io_surface.c
mex -output fastmarch1_mex -I./one2all one2all/my_heap.cpp one2all/unfold.cpp one2all/marchAlg.cpp one2all/io_surface.cpp
% mex -output fastmarch_mex  -I./all2all all2all/my_heap.c all2all/unfold.c all2all/marchAlg.c all2all/io_surface.c
% mex -output fastmarch_mex  -I./all2all all2all/my_heap.c all2all/marchAlg.c all2all/io_surface.c all2all/unfold.c

%%
clear all;
close all;
clc

load michael3;
surface=shape;
% One-to-all
source = repmat(Inf, size(surface.X));
source(1000)=0;

source(1)=0;
d = fastmarch(surface.TRIV, surface.X, surface.Y, surface.Z, source, set_options('mode', 'single'));
trisurf(surface.TRIV,surface.X,surface.Y,surface.Z,d(:,1)), axis equal, axis off
