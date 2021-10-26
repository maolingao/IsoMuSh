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

% wrapper to run ZoomOut for all pairs of a specific shape collection
%
function [C, C_map, time_zoomout] = zoomout_wrapper(data, params)

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
    S_map = data.S_map;
    
catch ME
    
    fprintf('** Fail to parse dataset! **\n');
    rethrow(ME)
end



numPair         = numShape * (numShape-1) / 2;
T_map           = cell(numPair, 2);  
C_map           = cell(numPair, 2); 
T               = cell(numPair, 2); 
C               = cell(numPair, 2);
time_zoomout    = nan(numPair, 2);


idx = 0;
params_zo.k_init    = 4;
params_zo.k_step    = 1;
params_zo.k_final   = dimLB;

%% ZoomOut loop

% for all source shapes
for i = 1 : numShape
    
    S1              = subsample_shape(S{i});
    s1_idx          = extract_number_from_string(S_map{i});
    B1              = S1.evecs(:, 1:params_zo.k_init);
    
    % for all target shapes 
    for j = (i+1) : numShape
        
        idx         = idx+1;
        S2          = subsample_shape(S{j});
        s2_idx      = extract_number_from_string(S_map{j});
        B2          = S2.evecs(:, 1:params_zo.k_init);
        
        %
        T_map{idx,1} = strcat('T_', s1_idx, '_', s2_idx);
        C_map{idx,1} = strcat('C_', s2_idx, '_', s1_idx);
        T_map{idx,2} = strcat('T_', s2_idx, '_', s1_idx);
        C_map{idx,2} = strcat('C_', s1_idx, '_', s2_idx);
        %
        % ZoomOut pair i->j
        t_start = tic;
        C21_ini = init_C(S{i}.origRes, S{j}.origRes, params_zo.k_init);
        T12_ini = fMAP.fMap2pMap(B2,B1,C21_ini);
        
        if isfield(S1,'sym')
            [T12, C21] = zoomOut_refine(S1.evecs, S2.evecs, T12_ini, params_zo, S1.sym, S2.sym, S1, S2);
        else
            [T12, C21] = zoomOut_refine(S1.evecs, S2.evecs, T12_ini, params_zo);
        end
        time_zoomout(idx, 1) = toc(t_start);
        
        if (verbose >= 1)
            fprintf('ZoomOut runtime %s -> %s : %.2f sec\n', S_map{i}, S_map{j}, time_zoomout(idx,1));
        end
        
        T{idx,1} = T12;
        C{idx,1} = C21;
        
        % ZoomOut pair j->i
        t_start = tic;
        C12_ini = init_C(S{j}.origRes, S{i}.origRes, params_zo.k_init);
        T21_ini = fMAP.fMap2pMap(B1,B2,C12_ini);
        
        if isfield(S1,'sym')
            [T21, C12] = zoomOut_refine(S2.evecs, S1.evecs, T21_ini, params_zo, S2.sym, S1.sym, S2, S1);
        else
            [T21, C12] = zoomOut_refine(S2.evecs, S1.evecs, T21_ini, params_zo);
        end
        
        time_zoomout(idx,2) = toc(t_start);
        
        if (verbose >= 1)
            fprintf('ZoomOut runtime %s -> %s : %.2f sec\n', S_map{j}, S_map{i}, time_zoomout(idx, 2));
        end
        
        T{idx,2} = T21;
        C{idx,2} = C12;
    end
end


fprintf('ZoomOut total runtime : %.2f sec\n', sum(time_zoomout, 'all'));

%% save zoomout results
datafullnameOut = strcat(dataname,'_zo.mat');
folder = fullfile(rootpath, 'data_init', datasetname);
if ~isfolder(folder)
    mkdir(folder)
end

save(fullfile(folder,datafullnameOut), 'T_map', 'C_map', 'T', 'C');
fprintf('Saved ZoomOut results (%s, %s) under %s.\n', datasetname, dataname, folder);

%% visualisation
% plot the final matching of the first 5 pairs
if (verbose >= 2)
    
    figure;
    subplot(2,3,1); 
    visualize_map_on_source(subsample_shape(S{2}), subsample_shape(S{1}), T{1}); 
    title(strcat('ZoomOut-', dataname, ' color legend'));
    
    for i = 2:6
        subplot(2,3,i); 
        visualize_map_on_target(subsample_shape(S{i}), subsample_shape(S{1}), T{i-1,2});
        if ( i >= numShape )
            break;
        end
    end

end
drawnow


end