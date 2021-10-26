% This function is part of the NmfSync algorithm for the synchronisation of
% partial permutation matrices, as described in [1]. When you use this
% code, you are required to cite [1].
% 
% [1] Synchronisation of Partial Multi-Matchings via Non-negative
% Factorisations F. Bernard, J. Thunberg, J. Goncalves, C. Theobalt.
% Pattern Recognition. 2019
%
% 
% Author & Copyright (C) 2019: Florian Bernard (f.bernardpi[at]gmail[dot]com)
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

function cycleError = cycleError(W, dimVector)
    Wcell = mat2cell(double(W),dimVector,dimVector);
    
    k = numel(dimVector);
    
    cycleError = 0;
    for l=1:k
        for i=1:k
            Til = Wcell{i,l};
            for j=1:k
                Tij = Wcell{i,j};
                Tlj = Wcell{l,j};

                availableI = find(sum(Til,2));
                availableJ = find(sum(Tlj,1));
                
                cycleError = cycleError + ...
                    norm(Til(availableI,:)*Tlj(:,availableJ)-...
                    Tij(availableI,availableJ), 'fro');
            end
        end
    end
    cycleError = cycleError/k^3;
end