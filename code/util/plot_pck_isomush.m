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
function plot_pck_isomush(datasetname, dataname, errorfolder, params)


    if ~exist(errorfolder,'dir')
        error('*** Error folder missing: %s\n', errorfolder);
    end
    
    folder = fileparts(errorfolder);
    
    verbose = params.verbose;
    
    if verbose >= 2
        fh = figure('Position', get(0, 'Screensize'),'visible','on');
    else
        fh = figure('Position', get(0, 'Screensize'),'visible','off');
    end
    
    thresh = 0:0.005:1;
%     thresh = 0:0.005:0.25;
    
    % load the pairs for pck evaluation
    switch datasetname

        case 'tosca'

            load('list_couples_tosca.mat','list_animal_name','list_num_1','list_num_2');
            ind_dataname = cellfun(@(x)strcmp(x, dataname), list_animal_name, 'UniformOutput', 1);
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
            idx_pair =  find(ind_dataname)';
        
        case 'scape'

            load('list_couples_scape.mat','couples_list_fm','list_pose_1','list_pose_2');
            num_pair = numel(couples_list_fm);
            idx_pair = 1:num_pair;
        

        otherwise
            error('Unknown datasetname');
    end

    %compute match errors
    curves = [];
    for k = idx_pair

        switch datasetname
            case 'tosca'
                id1 = list_num_1{k};
                id2 = list_num_2{k};
                name1 = strcat(list_animal_name{k}, id1);
                name2 = strcat(list_animal_name{k}, id2);
            case 'faust'
                name1 = list_pose_1{k};
                name2 = list_pose_2{k};
            case 'scape'
                name1 = list_pose_1{k};
                name2 = list_pose_2{k};
        end

        fname = strcat(name1,'_to_', name2);

        if exist(fullfile(errorfolder, strcat(fname, '.mat')), 'file')
            load(fullfile(errorfolder, strcat(fname, '.mat')), 'errors');
            if isempty(errors)
                continue;
            end

            % make plot
            num_error = min(100, size(errors,1));
            curve = calc_err_curve(errors(1:num_error), thresh);

            curves = [curves; curve];
        end
    end
        
    final_curve = mean(curves, 1);

    plot(thresh, final_curve, 'LineWidth', 4); ylim([50,100]);
    save(fullfile(folder,sprintf('%s_pck_isomush.mat', dataname)), 'thresh', 'final_curve');
    
    title(strcat(sprintf('pck-%s-%s',datasetname, dataname)), 'FontSize', 16);
    xlabel('Geodesic error', 'FontSize', 14); 
    ylabel('% Correspondences', 'FontSize', 14);
    
    saveas(gcf, fullfile(folder, sprintf('%s_pck_isomush.png', dataname)));
    
    if verbose==0
        close(fh);
    end
end