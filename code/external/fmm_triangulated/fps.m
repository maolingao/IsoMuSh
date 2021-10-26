function [sample,D]=fps(shape,N,init,sample_points)


D      = repmat(Inf, [length(shape.X) N]);    % Distance maps.
d      = repmat(Inf, [length(shape.X) 1]);

for k = 1:N
    if nargin<4
        if k==1
            sample = init;
        else
            sample = [sample; idx];
        end
    else
        sample = sample_points(k);
    end
    % Compute distance map from sample on the shape.
    u = repmat(Inf, [length(shape.X) 1]);
    u(sample(end)) = 0;
    D(:,k) = fastmarch(shape.TRIV, shape.X, shape.Y, shape.Z, u, set_options('mode', 'single'));
    
    %d = min(D,[],2);
    d = min(d, D(:,k));
    [r, idx] = max(d);
    
    % Visualize sampling in progress
    
end
