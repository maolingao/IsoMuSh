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

function P12 = T2perm(T12,nv2)
% T12 is the mapping from shape 1 to shape 2;
% P12 is the permutation matrix, which left-multiplies with a column vector
% of (1:nv_in_shape1)' and gives the T21.
nv1 = length(T12);
I = speye(nv2);
T12_zeros_removed = T12(T12~=0);
valid_map_idx = T12>0;

if length(unique(T12_zeros_removed)) ~= length(T12_zeros_removed)
%     fprintf('permError: many-to-one-mapping can NOT be converted into permutation matrix!\n');
%     fprintf('permError: convert to a loose permutation matrix instead!\n');
end

P21 = zeros(nv1,nv2);
P21(valid_map_idx,:) = I(T12_zeros_removed,:);
P12 = P21';
end