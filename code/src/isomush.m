% This function is part of the isoMuSh algorithm for the coorespondence 
% estimation of isometric shape collection, as described in [1]. When you
% use this code, you are required to cite [1].
% 
% [1] Isometric Multi-Shape Matching
% Author: M. Gao, Z. Lähner, J. Thunberg, D. Cremers, F. Bernard.
% IEEE Conference on Computer Vision and Pattern Recognition 2021 (CVPR 2021)
%
% 
% Author & Copyright (C) 2021: Maolin Gao (maolin.gao[at]gmail[dot]com)
%                              Florian Bernard (f.bernardpi[at]gmail[dot]com)
%
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.

% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

function [U, Q] = isomush(U, Q, data, params)

try
    rootpath = params.rootpath;
    verbose = params.verbose;
    datasetname = params.datasetname;
    dataname = params.dataname;
    maxIter = params.isomush.maxIter;
    
catch ME
    
    fprintf('** Fail to parse params! **\n');
    rethrow(ME)
end


try
    dimVector = data.dimVector;
    dimLB = data.dimLB;
    numShape = data.numShape;
    S = data.S;
    S_map = data.S_map;
    phicell = data.phicell;
    
catch ME
    
    fprintf('** Fail to parse dataset! **\n');
    rethrow(ME)
end

eps = 1e-6;
largeNoLAP = 1e3;
approxLap = 0;

%% IsoMuSh

if (verbose >= 1)
    fprintf('Starting IsoMuSh iteration...\n');
end
Phi = blkdiag(phicell{:});
obj_func = @(U, Q) trace(Q' * Phi' * U * U' * Phi * Q) / (numShape * numShape);

allObjs = obj_func(U, Q);
if (verbose >= 1)
    fprintf('iter: 0 obj: %.4e\n', allObjs);
end

t_start = tic;
for it = 1 : maxIter

    % U update: proj_p(Phi * Q * Q' * Phi' * U)
    PhiQ = Phi*Q;
    U2 = PhiQ*(PhiQ'*U);
    U = projectOntoPartialPermBlockwise(U2*largeNoLAP, dimVector, [], approxLap);

    % cost
    allObjs(2*it) = obj_func(U, Q);

    % Q update: proj_O(Phi' * U * U' * Phi' * Q)
    UTPhi = U'*Phi;
    Q2 = UTPhi'*(UTPhi*Q);
    Q = projectOntoStiefelBlockwise(Q2, numShape, 0);

    % cost
    allObjs(2*it+1) = obj_func(U, Q);
    if (verbose >= 1)
        fprintf('iter: %d obj: %.4e ', it, allObjs(2*it+1));
    end


    % stopping criteria
    relImprovement = 1 - (allObjs(2 * it) / allObjs(2 * it + 1));
    if (verbose >= 1)
        fprintf('relative improvement: %.4e\n', relImprovement);
    end


    if ( it > 1 && relImprovement - eps <= 0 )
        if (verbose >= 1)
            fprintf('Converged\n');
        end
        break;
    end

end


time = toc(t_start);
fprintf('Done:%.2f sec\n', time);
    

%% save isomush results
datafullnameOut = strcat(dataname,'_isomush.mat');
folder = fullfile(rootpath, 'result', datasetname);
if ~isfolder(folder)
    mkdir(folder)
end

save(fullfile(folder,datafullnameOut), 'U', 'Q');
fprintf('Saved IsoMuSh results (%s, %s) under %s.\n', datasetname, dataname, folder);

% save correspondences
[~, P] = decode_permutation_to_universe(U', dimVector);

matchfolder = fullfile(folder, 'matchings');
save_matches_isomush(P, dimVector, S, S_map, matchfolder);
fprintf('Saved IsoMuSh matches (%s, %s) under the path: %s.\n', datasetname, dataname, matchfolder);


%% visualisation
% plot the final matching of the first 5 pairs
if (verbose >= 2)
    
    figure;
    subplot(2,3,1); 
    visualize_map_on_source(subsample_shape(S{2}), subsample_shape(S{1}), perm2T(P{2,1})); 
    title(strcat('IsoMuSh-', dataname, ' color legend'));
    
    for i = 2:6
        subplot(2,3,i); 
        visualize_map_on_target(subsample_shape(S{i}), subsample_shape(S{1}), perm2T(P{i,1}));
        if ( i >= numShape )
            break;
        end
    end

end
drawnow

%% evaluate pck
% pck
errorfolder = fullfile(folder, 'errors');
evaluate_pck_isomush(datasetname, dataname, matchfolder, params)
plot_pck_isomush(datasetname, dataname, errorfolder, params)

end