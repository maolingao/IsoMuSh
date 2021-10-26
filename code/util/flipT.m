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
function T21 = flipT(T12,nv2)
% flip the mapping direction

P12 = T2perm(T12,nv2);

if (sum(unique(sum(P12,2))>=2) >=1 ) % row-sum >1, shape2->shape1 contains many-to-one mapping
%     fprintf('permError: Not a permutation matrix!\n');
%     fprintf('permError: randomly pick one correspondence!\n');
    % random sample a match
    T21 = nan(nv2,1);
    for i=1:nv2
        matches = find(P12(i,:));
        if isempty(matches)
            T21(i) = 0;
        else
            T21(i) = matches(randi(length(matches)));
        end
    end
else % perfect permutation case
    idxVec = (1:1:size(P12,2))';
    T21 = P12*idxVec;
end


end