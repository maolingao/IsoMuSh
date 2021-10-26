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
%                              Zorah Lähner (zorah.laehner[at]uni-siegen[dot]de)
%
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
function evaluate_pck_isomush(datasetname, dataname, matchfolder, params)


if ~exist(matchfolder, 'dir')
    fprintf('*** Matching folder missing: %s\n', matchfolder);
    return
end

% folder to save the geodesic error
errorfolder = fileparts(matchfolder);
errorfolder = fullfile(errorfolder, 'errors');
if ~exist(errorfolder, 'dir')
    mkdir(errorfolder);
end

datapath_orig = params.dataPreparation.datapath_orig;


% load the pairs for pck evaluation
switch datasetname
    
    case 'tosca'
        
        load('list_couples_tosca.mat','list_animal_name','list_num_1','list_num_2');
        ind_dataname = cellfun(@(x)strcmp(x, dataname), list_animal_name, 'UniformOutput', 1);
        num_pair = sum(ind_dataname);
        idx_pair =  find(ind_dataname)';
        
    case 'faust'
        
        load('list_couples_faust.mat','couples_list_fm','list_pose_1','list_pose_2');
        % cell array of indices of the same person
        idx_pose = num2cell(0:9);
        idx_pose = cellfun(@(x)(num2str(x)), idx_pose, 'UniformOutput', 0);
        idx_pose = string(cellfun(@(x)(strcat('0', dataname, x)), idx_pose, 'UniformOutput', 0));
        % pairs belong to the same person
        ind_dataname_pose1 = cellfun(@(x) any(contains(idx_pose, extract_number_from_string(x))),...
            list_pose_1, 'UniformOutput', 1);
        ind_dataname_pose2 = cellfun(@(x) any(contains(idx_pose, extract_number_from_string(x))),...
            list_pose_2, 'UniformOutput', 1);
        ind_dataname = ind_dataname_pose1 & ind_dataname_pose2;
        %
        num_pair = numel(ind_dataname);
        idx_pair =  find(ind_dataname)';
        
    case 'scape'
        
        load('list_couples_scape.mat','couples_list_fm','list_pose_1','list_pose_2');
        num_pair = numel(couples_list_fm);
        idx_pair = 1:num_pair;
        
    otherwise
        error('Unknown datasetname');
end

fprintf('Loaded pairs to be evaluated from: list_couples_%s\n', datasetname);



n_samples = 1000;
kk = 1;
%compute match errors
for k = idx_pair
    fprintf('-------------------------------------\nProcessing shape %d of %d\n', kk, num_pair);
    clear X;
    clear Y;
    
    switch datasetname
        case {'faust', 'scape'}
            
            name1 = list_pose_1{k};
            name2 = list_pose_2{k};
            fname = strcat(name1,'_to_', name2);

        case 'tosca'
            
            id = k;
            id1 = list_num_1{id};
            id2 = list_num_2{id};
            name1 = strcat(list_animal_name{id}, id1);
            name2 = strcat(list_animal_name{id}, id2);
            fname = strcat(list_animal_name{id}, id1,'_to_', list_animal_name{id}, id2);
    end
    
    if ~exist(fullfile(matchfolder, strcat(fname, '.mat')), 'file')
        fprintf('Skipping %s.\n', fname);
        continue;
    end

    switch datasetname
        case 'faust'
            
            X = load_ply(fullfile(datapath_orig, strcat(name1, '.ply')));
            Y = load_ply(fullfile(datapath_orig, strcat(name2, '.ply')));
        
        case {'scape', 'tosca'}
            
            X = load_off(fullfile(datapath_orig, strcat(name1, '.off')));
            Y = load_off(fullfile(datapath_orig, strcat(name2, '.off')));
        
    end
    
    X.n = size(X.VERT, 1);
    X.m = size(X.TRIV, 2);
    Y.n = size(Y.VERT, 1);
    Y.m = size(Y.TRIV, 2);


    fprintf('Evaluating %s %s.\n', name1, name2);

    % load matches
    load(fullfile(matchfolder, strcat(fname, '.mat')),'matches');

    % calculate errors
    fprintf("Calculating errors... \n");

    n_samples_here = min(n_samples, size(matches,1));
    if n_samples_here < length(matches)
        % choose random samples
        samples = randsample(length(matches), n_samples_here);
    else
        samples = 1:length(matches);
    end
    errors = zeros(length(samples), 1);
    for i=1:length(samples)
      if matches(i,1) == matches(i,2)
           continue;
       end
       src = Inf*ones(Y.n, 1);
       src(matches(i,1)) = 0;
       d = fast_marching(Y, src);
       errors(i) = d(matches(i,2));
    end
    % Farthest points extraction
    fps_samples = 50;
    Y.fps = fps_euclidean(Y.VERT, fps_samples, 1);
    D = fastmarch_multi(Y, Y.fps);
    diameter = max(D,[],'all');
    
    errors = errors ./ diameter;

    save(fullfile(errorfolder, strcat(fname, '.mat')), 'errors');

    fprintf('Mean error %.3f\n',mean(errors));
    kk = kk + 1;
end
fprintf('saved errors of %s under the path: %s.\n', datasetname, errorfolder);

end