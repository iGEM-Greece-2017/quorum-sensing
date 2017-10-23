function K=checkForSymmetry(K)
%checkForSymmetry - Check symmetry for K matrix

% Copyright 2016 The MathWorks, Inc

% If the K matrix is very close to symmetric, make it exactly symmetric
% so \ will use the symmetric factorization routine.
if ~isempty(K)
    maxDiagK = full(max(abs(diag(K))));
    maxAs = full(max(max(abs(K-K'))));
    if(maxAs < eps(maxDiagK))
        K = (K+K')/2;
    end
end
end