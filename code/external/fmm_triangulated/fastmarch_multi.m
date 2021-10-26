% this is a function to compute the geodesic distances btw a subset of
% vertices on the triangulated manifold M
% developed by Maolin Gao 
function D = fastmarch_multi(M, samples)
ns = length(samples);
D = nan(ns,ns);

for i = 1:ns
   src = Inf*ones(M.n, 1);
   src(samples(i)) = 0;
   d_full = fast_marching(M, src);
   d = d_full(samples);
   D(i,:) = d';
end

end