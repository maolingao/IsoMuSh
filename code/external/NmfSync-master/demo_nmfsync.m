% load simple demo instance
load demo_data

% run NmfSync
d = size(U_gt,2);
[Wsync, U] = nmfSync(Wnoisy, dimVector, d);

% functions that compute errors
tp = @(W_gt,W) nnz(W_gt & W);
fp = @(W_gt,W) nnz(~W_gt & W);
fn = @(W_gt,W) nnz(W_gt & ~W);
prec = @(W_gt,W) tp(W_gt,W)/(tp(W_gt,W) + fp(W_gt,W));
rec = @(W_gt,W) tp(W_gt,W)/(tp(W_gt,W) + fn(W_gt,W));
fscore = @(W_gt,W) (2*prec(W_gt,W)*rec(W_gt,W))/(prec(W_gt,W)+rec(W_gt,W));

% evaluate pairwise matchings
cycleErrorNoisy = cycleError(Wnoisy, dimVector);
fscoreNoisy = fscore(W_gt,Wnoisy);

cycleErrorSync = cycleError(Wsync, dimVector);
fscoreSync = fscore(W_gt,Wsync);

disp(['Before NmfSync: Cycle-error = ' num2str(cycleErrorNoisy) ', fscore = ' num2str(fscoreNoisy)]);
disp(['After  NmfSync: Cycle-error = ' num2str(cycleErrorSync) ', fscore = ' num2str(fscoreSync)]);


