function [T12, C21, all_T12, all_C21] = zoomOut_refine_partial(B1_all, B2_all, T, para, sym1, sym2, rect)

if nargout > 2, all_T12 = {}; all_C21 = {}; end

for k = para.k_init : para.k_step : para.k_final
    B1 = B1_all(:, 1:k+rect);
    B2 = B2_all(:, 1:k);
    C21 = B1(T,:)\B2;
    if nargin==4
        % original zoomout
        T12 = knnsearch(B2*C21', B1);
        
    else
        % add symmetric descriptors into B1,B2
%         sym11= sym1 ./ norm(sym1,'Fro'); % TODO: if you want to be super clean sym1 ./ (sym1'M1sym1), M1 is the mass matrix,
%         sym22= sym2 ./ norm(sym2,'Fro');
        T12 = knnsearch([B2*C21',sym2], [B1,sym1]);
    end
    
    if nargout > 2
        all_T12{end+1} = T12; all_C21{end+1} = C21;
    end

end