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

function save_matches_isomush(X, dimVector, S, S_map, matchfolder)
% X could be 
% 1. Prel = a cell array with rows = columns = num of shapes in the
% collection, and each cell contains Pij, which satisfies Tij=Pij*(1:nvj)'
% 2. T = a cell array of Tij (i<j only).
% 3. Tmatrix = a cell array with row = columns = num of shapes in the
% collection, and each cell contains Tij.
% dimVector is the vector of all m_i, i = {1,...,k}
% S_map is the dataname of each pose


if ~isfolder(matchfolder)
    mkdir(matchfolder);
end

k = length(dimVector);

if size(X, 1) == size(X, 2)
    if size(X{2}, 1)==1 || size(X{2}, 2)==1 % Tmatrix
        Tmatrix = X;
    else % Prel
        Tmatrix = convert_Prel_to_Tmatrix(X);
    end
else % T (either half or full)
    Tmatrix = convert_T_to_Tmatrix(X, dimVector);
end


% save matches_subsampled
for i = 1 : k
    for j = (i+1) : k
        %
        matches_subsampled = nan(dimVector(i),2);
        orig_idx  = cell(2,1);
        filename = strcat(S_map{i},'_to_',S_map{j});
        matches_subsampled(:,1) = (1:dimVector(i))';
%         matches_subsampled(:,2) = perm2T(Prel{i,j});
        matches_subsampled(:,2) = Tmatrix{i,j};
        orig_idx{1} = subsample_shape(S{i}).orig_idx;
        orig_idx{2} = subsample_shape(S{j}).orig_idx;
        matches = convert_matches_to_original_res(matches_subsampled, orig_idx);
        save(strcat(fullfile(matchfolder,filename),'.mat'), 'matches','matches_subsampled','orig_idx');
        %
        matches_subsampled = nan(dimVector(j),2);
        orig_idx  = cell(2,1);
        filename = strcat(S_map{j},'_to_',S_map{i});
        matches_subsampled(:,1) = (1:dimVector(j))';
%         matches_subsampled(:,2) = perm2T(Prel{j,i});
        matches_subsampled(:,2) = Tmatrix{j,i};
        orig_idx{1} = subsample_shape(S{j}).orig_idx;
        orig_idx{2} = subsample_shape(S{i}).orig_idx;
        matches = convert_matches_to_original_res(matches_subsampled, orig_idx);
        save(strcat(fullfile(matchfolder,filename),'.mat'), 'matches','matches_subsampled','orig_idx');
    end
end

end


function Tmatrix = convert_T_to_Tmatrix(T,dimVector)
k = length(dimVector);
Tmatrix = cell(k,k);
idx = 1;
for i=1:k
    for j=(i+1):k
        
    Tmatrix(i,j) = T(idx);
    if size(T,1)==1 || size(T,2)==1
        % only half of pairs available as a vector
        Tmatrix(j,i) = {flipT(T{idx,1},dimVector(j))};
    else
        % all pairs available as 2 vectors
        Tmatrix(j,i) = T(idx,2);
    end
    idx = idx+1;
    
    end
end

        
end