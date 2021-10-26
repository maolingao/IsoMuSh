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

function [S, S_map] = load_dataset_origFormat(params)

try
    rootpath = params.rootpath;
    verbose = params.verbose;
    datasetname = params.datasetname;
    dataname = params.dataname;
    dimLB = params.dataPreparation.dimLB;
    numFace = params.dataPreparation.numFace;
    datapath_orig = params.dataPreparation.datapath_orig;
    datapath_mat = params.dataPreparation.datapath_mat;
    
catch ME
    
    fprintf('** Fail to parse params! **\n');
    rethrow(ME)
end

% output path for preprocessed mat files
if ~isfolder(datapath_mat)
    mkdir(datapath_mat)
end

try
    % load from mat in original resolution, subsample faces and randomly
    % shuffle vertices
    [S, S_map] = subsample([], [], numFace, params);

catch
    % mat in original resolution are not available
    % load from the original *.off, *.ply data
    switch datasetname
        case 'tosca'
            symdata = fullfile(rootpath, 'data/tosca/sym', ...
                strcat(dataname, '.sym.labels'));
            sym = load_sym(symdata);


            filelist = dir(datapath_orig);
            filenames = {filelist(:).name};
            shapeidx = cellfun(@(x)contains(x, dataname), filenames, 'UniformOutput', 1);

            shapenames = filenames(shapeidx);
            numShape = length(shapenames);
            S = cell(numShape,1);
            S_map = cell(numShape,1);

            for j = 1 : numShape

                [~,shapename,~] = fileparts(shapenames{j});
                S_map{j} = shapename;

                tmp = load_off(fullfile(datapath_orig, shapenames{j}));
                S{j} = convert_struct(tmp);
                S{j}.name = shapename;
                S{j}.sym = sym;

                S{j} = MESH.compute_LaplacianBasis(S{j}, dimLB);

            end

            save(strcat(datapath_mat,'/',dataname,'_origRes.mat'), 'S', 'S_map', '-v7.3');




        case 'faust'
            symdata = fullfile(rootpath, 'data/faust/sym/faust_registrations.sym.labels');
            sym = load_sym(symdata);

            S = cell(10,1); % each person in Faust has 10 poses.
            S_map = cell(10,1);
            
            for j = 0 : 9
                    shapename = strcat('tr_reg_0', dataname, num2str(j));
                    datafullpath = fullfile(datapath_orig, strcat(shapename, '.ply'));
                    M = load_ply(datafullpath);
                    S_map{j+1} = shapename;
                    Sj = convert_struct(M);
                    Sj.name = shapename;
                    Sj.sym = sym;
                    %
                    Sj = MESH.compute_LaplacianBasis(Sj, dimLB);
                    S{j+1} = Sj;
                    %
                    % figure, trisurf(TRI,PTS(:,1),PTS(:,2),PTS(:,3)); axis off equal;
                    % figure, visualize_map_on_source('', Sj, '', Sj.sym);
            end
            
            save(fullfile(datapath_mat, strcat(dataname, '_origRes.mat')), 'S', 'S_map');



        case 'scape'

            symdata = fullfile(rootpath, 'data/scape/sym/scape.sym.labels');
            sym = load_sym(symdata);
            load(fullfile(rootpath, 'eval/list_couples_scape.mat'), 'sampled_mesh_names');
            numShape = length(sampled_mesh_names);
            idx = 1;
            S = cell(numShape, 1);
            S_map = cell(numShape, 1);

            for j = 0 : 71

                if j <= 9
                    shapename = strcat('mesh0', '0', num2str(j));
                else
                    shapename = strcat('mesh0', num2str(j));
                end
                if ~sum(strcmp(sampled_mesh_names, shapename))
                    continue
                end

                datafullpath = fullfile(datapath_orig, strcat(shapename, '.off'));
                M = load_off(datafullpath);
                S_map{idx} = shapename;
                Sj = convert_struct(M);
                Sj.name = shapename;
                Sj.sym = sym;
                %
                Sj = MESH.compute_LaplacianBasis(Sj, dimLB);
                S{idx} = Sj;
                idx = idx+1;
                %
                % figure, trisurf(TRI,PTS(:,1),PTS(:,2),PTS(:,3)); axis off equal;
                % figure, visualize_map_on_source('', Sj, '', Sj.sym);
            end

            save(fullfile(datapath_mat, strcat(dataname, '_origRes.mat')), 'S', 'S_map');


        otherwise
            error('unsupported datasetname!')
    end


    [S] = subsample(S, S_map, numFace, params);
end
end


function S = convert_struct(M)

S.surface.TRIV = M.TRIV;
S.surface.VERT = M.VERT;
S.surface.X = S.surface.VERT(:,1);
S.surface.Y = S.surface.VERT(:,2);
S.surface.Z = S.surface.VERT(:,3);

S.nv = size(S.surface.VERT,1);
S.nf = size(S.surface.TRIV,1);

end
