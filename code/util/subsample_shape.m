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
function [S_output] = subsample_shape(S)

if ~isfield(S, 'origRes')
    % data is already full resolution
    S_output = S;
    if ~isfield(S,'surface') || ~isfield(S.surface,'X')
        S_output.surface.TRIV = S_output.TRIV;
        S_output.surface.VERT = S_output.VERT;
        S_output.surface.X = S_output.VERT(:,1);
        S_output.surface.Y = S_output.VERT(:,2);
        S_output.surface.Z = S_output.VERT(:,3);
    end
else % subsample data
    % subsampling using the orig_idx and faces
    orig_idx = S.orig_idx;
    S_subsampled.orig_idx = orig_idx;
    S_subsampled.surface.TRIV = S.TRIV_subsampled;
    S_subsampled.surface.VERT = S.origRes.surface.VERT(orig_idx,:);
    S_subsampled.surface.X = S_subsampled.surface.VERT(:,1);
    S_subsampled.surface.Y = S_subsampled.surface.VERT(:,2);
    S_subsampled.surface.Z = S_subsampled.surface.VERT(:,3);
    S_subsampled.nv = size(S_subsampled.surface.VERT,1);
    S_subsampled.nf = size(S_subsampled.surface.TRIV,1);
    S_subsampled.evecs = S.origRes.evecs(orig_idx,:);
    if isfield(S.origRes,'sym')
        S_subsampled.sym = S.origRes.sym(orig_idx);
    end
    S_subsampled.wks = S.origRes.wks(orig_idx,:);
    S_subsampled.hks = S.origRes.hks(orig_idx,:);
    S_subsampled.desc = S.origRes.desc(orig_idx,:);
    %
    S_subsampled.evals = S.origRes.evals;
    S_subsampled.name = S.origRes.name;
    S_subsampled.shuffle_maps = S.shuffle_maps;
    S_subsampled.A = S.origRes.A(orig_idx,orig_idx);
    
    shuffle=1;
    if shuffle==1
        shuffle_maps = S_subsampled.shuffle_maps;
        S_shuffled = shuffle_vertices(S_subsampled, shuffle_maps);
        S_output = S_shuffled;
    else
        S_output = S_subsampled;
    end
%     S_subsampled.W = S.origRes.W;
%     S_subsampled.area = S.origRes.area;
%     S_subsampled.sqrt_area = S.origRes.sqrt_area;
end

end

function S_shuffled = shuffle_vertices(S,shuffle_maps)
T21 = shuffle_maps(:,2); % map what is the new vertex in position 1
T12 = shuffle_maps(:,3); % where is the previous 1st vertex now
S_shuffled.orig_idx = S.orig_idx(T21);
S_shuffled.surface.TRIV = T12(S.surface.TRIV);
S_shuffled.surface.VERT = S.surface.VERT(T21,:);
S_shuffled.surface.X = S.surface.X(T21);
S_shuffled.surface.Y = S.surface.Y(T21);
S_shuffled.surface.Z = S.surface.Z(T21);
S_shuffled.evecs = S.evecs(T21,:);
if isfield(S,'sym')
    S_shuffled.sym = S.sym(T21);
end
S_shuffled.wks = S.wks(T21,:);
S_shuffled.hks = S.hks(T21,:);
S_shuffled.desc = S.desc(T21,:);
S_shuffled.A = S.A(T21,T21);

%
S_shuffled.nv = S.nv;
S_shuffled.nf = S.nf;
S_shuffled.evals = S.evals;
S_shuffled.name = S.name;

end