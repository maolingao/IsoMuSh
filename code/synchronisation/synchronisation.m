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

function [U, Q] = synchronisation(C, data, params)


try
    rootpath = params.rootpath;
    verbose = params.verbose;
    datasetname = params.datasetname;
    dataname = params.dataname;
    
catch ME
    
    fprintf('** Fail to parse params! **\n');
    rethrow(ME)
end


try
    dimVector = data.dimVector;
    dimLB = data.dimLB;
    numShape = data.numShape;
    S = data.S;
    
catch ME
    
    fprintf('** Fail to parse dataset! **\n');
    rethrow(ME)
end


%% Q synchronisation

Qcell   = cell(numShape);
I       = speye(dimLB);


idx=1;
for i = 1 : numShape
    
    Qcell{i,i} = I;
    
    for j = i+1 : numShape
        
        % Q_{ij}
        Cij         = C{idx};
        
        % band filter
        bandwidth   = 6;
        bands       = -bandwidth : bandwidth;
        M           = spdiags(ones(dimLB, numel(bands)), bands, dimLB, dimLB);
        Cij         = full(Cij.*M);
        
        Qcell{i,j}  = projectOntoStiefelBlockwise(Cij, 1, 0);
        
        % Q_{ji}
        Qcell{j,i}  = (Qcell{i,j})';
        
        idx = idx + 1;
    end
end

fprintf('Synchronising functional maps...');

% orthogonal transformation synchronisation
tstart = tic;
Q = syncQfromPairwise(Qcell);
time = toc(tstart);

fprintf('done:%.2f sec\n', time);


%% U synchronisition

Psi = [];
for i = 1 : numShape
    phii            = subsample_shape(S{i}).evecs;
    Qi              = Q(dimLB*(i-1)+1 : dimLB*i, :);
    phii_aligned    = phii * Qi;
    Psi             = [Psi; phii_aligned]; %#ok<AGROW>
end

fprintf('Synchronising permutation matrices...');

% Successive Block Rotation Algorithm
tstart = tic;
universeSize = max(dimVector);
Urot = SBRA(Psi*Psi', universeSize, dimVector);
Urot = Urot(:,1:universeSize);

largeNoLAP = 1e3;
approxLap = 0;
U = projectOntoPartialPermBlockwise(Urot'*largeNoLAP, [], dimVector, approxLap);
U = U';
time = toc(tstart);

fprintf('done:%.2f sec\n', time);


%% save synchronisation results
datafullnameOut = strcat(dataname,'_sync.mat');
folder = fullfile(rootpath, 'data_init', datasetname);
if ~isfolder(folder)
    mkdir(folder)
end

save(fullfile(folder,datafullnameOut), 'U', 'Q');
fprintf('Saved synchronisation results (%s, %s) under %s.\n', datasetname, dataname, folder);

%% visualisation
% plot the final matching of the first 5 pairs
if (verbose >= 2)
    
    [~, P] = decode_permutation_to_universe(U', dimVector);
    
    figure;
    subplot(2,3,1); 
    visualize_map_on_source(subsample_shape(S{2}), subsample_shape(S{1}), perm2T(P{2,1})); 
    title(strcat('Synch-', dataname, ' color legend'));
    
    for i = 2:6
        subplot(2,3,i); 
        visualize_map_on_target(subsample_shape(S{i}), subsample_shape(S{1}), perm2T(P{i,1}));
        if ( i >= numShape )
            break;
        end
    end

end
drawnow

end