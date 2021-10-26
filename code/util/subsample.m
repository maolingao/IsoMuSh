function [S, S_map] = subsample(S, S_map, numFace, params)
% note that only the script has been tested in the fowllowing datasets: 
% tosca, faust, scape as stated in the paper (except tosca wolf).
% 

try
    rootpath = params.rootpath;
    verbose = params.verbose;
    datasetname = params.datasetname;
    dataname = params.dataname;
    datapath_mat = params.dataPreparation.datapath_mat;
    
catch ME
    
    fprintf('** Fail to parse params! **\n');
    rethrow(ME)
end


%% core
% load shapes in original resolution with cached LB
if isempty(S) && isempty(S_map) % only load if S and S_map are not inputs
    load(strcat(datapath_mat, '/', dataname, '_origRes.mat'), 'S', 'S_map')
end
numShape = length(S);
S_origRes = S;
clear S;
S = cell(numShape,1);

% subsample the null shape and remember the sampled vertices
for i=1:numShape
    Si_origRes = S_origRes{i};
    idxstr = extract_number_from_string(Si_origRes.name);
    
    if (strcmp(idxstr, repmat('0', 1, length(idxstr)))) ||...
             (strcmp(dataname, 'gorilla') && length(idxstr)==1 && strcmp(idxstr(end),'1')) ||... % gorilla0 does not exist
             (strcmp(datasetname, 'faust') && length(idxstr)==3 && strcmp(idxstr(end),'0')) % faust shapes
         
        fprintf('Subsampling pose %s...\n', Si_origRes.name);
        [Si_origRes.shot, Si_origRes.hks, Si_origRes.desc, Si_origRes.wks] = ...
            cal_shot_hks_wks(Si_origRes);
        
        if (Si_origRes.nf > numFace)
            
            % save the full resolution in a field
            Si.origRes = Si_origRes;
            Si.name = Si_origRes.name;
            
            % subsampling
            fv.vertices = Si_origRes.surface.VERT;
            fv.faces = Si_origRes.surface.TRIV;
            fv_subsampled = reducepatch(fv, numFace);
            
            % save the subsampled vertices (orig_idx) and faces (TRIV_subsampled)
            subsampled_TRIV = fv_subsampled.faces; % the subsampled faces are the same across all shapes
            subsampled_vertices_idx = find_orig_idx(fv.vertices, fv_subsampled.vertices); % and the same subset of vertices
            Si.orig_idx = subsampled_vertices_idx;
            Si.TRIV_subsampled = subsampled_TRIV;
            
            % shuffle maps [position, T21, T12]
            nv = length(Si.orig_idx);
            Si.shuffle_maps(:,1) = (1:nv)';
            Si.shuffle_maps(:,2) = randperm(nv)';
            Si.shuffle_maps(:,3) = flipT(Si.shuffle_maps(:,2), nv);
        else
            
            Si = Si_origRes;
        end

        S{i} = Si;
        idx_null = i;
        clear Si_origRes;
    end

end


% subsample the same set of vertices as the null shape above
for i = 1 : numShape
    if i~=idx_null
        
        Si_origRes = S_origRes{i};
        fprintf('Subsampling %s...\n', Si_origRes.name);
        [Si_origRes.shot, Si_origRes.hks, Si_origRes.desc, Si_origRes.wks] = ...
            cal_shot_hks_wks(Si_origRes);
        
        if (Si_origRes.nf > numFace)
            
            % save the full resolution in a field
            Si.origRes = Si_origRes;
            Si.name = Si_origRes.name;
            
            % save same vertices and faces as shape null
            Si.orig_idx = subsampled_vertices_idx;
            Si.TRIV_subsampled = subsampled_TRIV;
            
            % shuffle maps [position, T21, T12]
            nv = length(Si.orig_idx);
            Si.shuffle_maps(:,1) = (1:nv)';
            Si.shuffle_maps(:,2) = randperm(nv)';
            Si.shuffle_maps(:,3) = flipT(Si.shuffle_maps(:,2),nv);
        else
            
            Si = Si_origRes;
        end
        
        S{i} = Si;
        clear Si_origRes;
    end
end


%% save
outputfolder = fullfile(rootpath, 'data', datasetname);
if ~isfolder(outputfolder)
   mkdir(outputfolder);
end

save(strcat(outputfolder, '/', dataname, '.mat'), 'S', 'S_map', '-v7.3');
fprintf('saved the subsampled %s mat file under the path: %s/.\n', dataname, outputfolder);

end