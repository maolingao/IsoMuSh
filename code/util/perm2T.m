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

function T12 = perm2T(P12)
% T12 is the mapping from shape 1 to shape 2;
% P12 is the permutation matrix, which left-multiplies with a column vector
% of (1:nv_in_shape1)' and gives the T21.
nv1 = size(P12,2);
nv2 = size(P12,1);
if (sum(unique(sum(P12,1))>=2) >=1 || sum(unique(sum(P12,2))>=2) >=1 ) % col-sum or row-sum >1
%     fprintf('permError: Not a permutation matrix!\n');
%     fprintf('permError: randomly pick one correspondence!\n');
    if (sum(unique(sum(P12,1))>=2) ==0) % col-sum <=1, shape1->shape2 only contains injective or no mapping
        idxVec = (1:nv2)';
        T12 = P12'*idxVec;
    else % shape1->shape2 and shape2->shape1 constains many-to-one mapping
        % random sample a match
        T12 = nan(nv1,1);
        for i=1:nv1
            matches = find(P12(:,i));
            if isempty(matches)
                T12(i) = 0;
            else
                T12(i) = matches(randi(length(matches)));
            end
        end
    end
else % perfect permutation case
    idxVec = (1:nv1)';
    T21 = P12*idxVec;
    T12 = flipT(T21,nv1);
end

end