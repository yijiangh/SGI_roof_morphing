function W = compute_skinning_weight(V,F,b)
%%%%%%%%%%
% TASK: implement (unbounded) biharmonic weights, which solves the
% following optimization problem: 
%      min   0.5 * trace( W'*Q*W )
%      s.t.  W(b,:) = bc
%
% black list:
% - min_quad_with_fixed
%
% white list (hints)
% - "cotmatrix" for computing cotangent Laplacian
% - "massmatrix" for computing mass matrix
%
% Hint:
% - give a diagonal matrix A, you can invert A with 
%   >> invA = diag(diag(A).^(-1));
% - constraints are 1 for a handle and 0 for other handls, so it is like
%   >> bc = eye(length(b));
% - you shouldn't need any other functions in gptoolbox
%%%%%%%%%%

% naive closest point (replace this one with your solution)
L = -cotmatrix(V,F);
M = massmatrix(V,F);

Q = L * diag(diag(M).^(-1)) * L;
W = zeros(length(V), length(b));
bc = eye(length(b));

not_b = ones(length(V), 1);
not_b(b) = false;
not_b = find(not_b);
Q11 = Q(not_b, not_b);
Q12 = Q(not_b, b);

for j=1:length(b)
    W(not_b,j) = -1 * Q11 \ (Q12 * bc(:,j));
    W(b, j) = bc(:,j);
end
    
    
end