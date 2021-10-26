function [X_shot,X_hks,X_desc,X_wks] = cal_shot_hks_wks(X, opts)

if nargin < 2
    opts = struct;
end
if ~isfield(opts, 'shot_num_bins')
    opts.shot_num_bins = 10; % number of bins for shot
end
if ~isfield(opts, 'shot_radius')
    opts.shot_radius = 5; % percentage of the diameter used for shot
end

X.n = X.nv;

X_shot = calc_shot(X.surface.VERT', X.surface.TRIV', 1:X.nv, opts.shot_num_bins, opts.shot_radius*sqrt(sum(X.area(:)))/100, 3)';
X_hks = calc_HKS(X,size(X.evecs,2));

X_desc = [X_shot./max(X_shot(:)), X_hks./max(X_hks(:))];

X_wks = waveKernelSignature(X.evecs, X.evals, X.A, 200);

end