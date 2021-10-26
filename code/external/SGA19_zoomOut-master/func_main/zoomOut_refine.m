% modified version of original ZoomOut
% additionally taking symetric descriptor into account

function [T12, C21, all_T12, all_C21] = zoomOut_refine(B1_all, B2_all, T12, para, sym1, sym2, S1, S2)

if nargout > 2, all_T12 = {}; all_C21 = {}; end

if nargin==8
    S1 = MESH.compute_LaplacianBasis(S1, 100);
    S2 = MESH.compute_LaplacianBasis(S2, 100);
    weight = 1e2;
    sym1= weight*sym1 ./ sqrt(sym1'*S1.A*sym1);
    sym2= weight*sym2 ./ sqrt(sym2'*S2.A*sym2);
end

for k = para.k_init : para.k_step : para.k_final
    B1 = B1_all(:, 1:k);
    B2 = B2_all(:, 1:k);
    C21 = B1\B2(T12,:);
    if nargin==4
        % original zoomout
        T12 = knnsearch(B2*C21', B1);
        
    elseif nargin==8
        % add symmetric descriptors into B1,B2
        T12 = knnsearch([B2*C21',sym2], [B1,sym1]);
    end
    
    if nargout > 2
        all_T12{end+1} = T12; all_C21{end+1} = C21;
    end

end