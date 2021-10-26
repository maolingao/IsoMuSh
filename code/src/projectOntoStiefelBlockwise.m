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


function R = projectOntoStiefelBlockwise(R2, k, detEqualsOne)
    D = size(R2,1)/k;
    
    R = [];
    for i=1:k
        idx = ((i-1)*D+1):(i*D);
        %         Ri = R(idx,:);
        R2i = R2(idx,:);
        [u,~,v] = svd(R2i, 'econ');
        if ( detEqualsOne )
            R(idx,:) = u*diag([ones(1, D-1), det(u*v')])*v';
        else
            R(idx,:) = u*v';
        end
    end
end