% This file implements the method described in:
%
% "Consistent Partial Matching of Shape Collections via Sparse Modeling" 
% L. Cosmo, E. Rodola, A. Albarelli, F. Memoli, D. Cremers
% Computer Graphics Forum, 2016
%
% If you use this code, please cite the paper above.
% 
% Luca Cosmo, Emanuele Rodola (c) 2015

function [ S ] = fps_euclidean(V, n, seed)
%% FPS_EUCLIDEAN Samples K vertices from V by using farthest point sampling.
% The farthest point sampling starts with vertex v1 and uses the euclidean
% metric of the 3-dimensional embedding space.
% -  V is a n-by-3 matrix storing the positions of n vertices
% -  K is the number of samples
% -  v1 is the index of the first vertex in the sample set. (1<=v1<=n)
% Returns
% -  S is a K-dimensional vector that includes the indeces of the K sample
%    vertices.

%Hint: matlab function pdist2 could be helpful

S = zeros(n,1);
S(1) = seed;
d = pdist2(V,V(seed,:));

for i=2:n
    [~,m] = max(d);
    S(i) = m(1);
    d = min(pdist2(V,V(S(i),:)),d);
end

end