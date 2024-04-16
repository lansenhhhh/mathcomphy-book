function x = mkconstarray(A, x)
    % Ensure x is a column vector of the correct size
    if numel(x) ~= A.n * A.m
        error('Size of vector x does not match expected dimensions.');
    end
    
    % Your existing operations...
    x = reshape(x, A.n, A.m); % Now safe to reshape
    % Followed by your operator's specific actions on x
end