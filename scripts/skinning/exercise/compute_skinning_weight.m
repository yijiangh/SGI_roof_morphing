function W = compute_skinning_weight(V,F,b)
%%%%%%%%%%
% Input:
%   V: #V x dim matrix of the original mesh vertex positions
%   F: #F x 3 matrix of original triangular mesh's connectivity
%   b: a list of handle vertex indices
% 
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
idx = knnsearch(V(b,:),V);
% init weight matrix W: #V x #handles
W = full(sparse(1:size(V,1), idx, 1, size(V,1), length(b)));

L = -cotmatrix(V, F);
M = massmatrix(V, F);
%  #V x #V bilaplacian operator
Q = L * (M\L);

% constraint W(b(j),j) = 1 on handles
% So for each handle, we are solving the quadratic problem:
% minimize    W(:,j)' * Q * W(:,j)
% subject to  W(b,j) = bc(:,j)

% #handles x #handles matrix of the constraints on the skinning weights
bc = eye(length(b));

% argmin_u min(u^T * Q * u) => u = [x; c], where x = Q(a,a) \ (-Q(a,b)*c)
% here a is the non-fix indices, b is the fixed indices
notb = ones(size(V,1),1);
notb(b) = false;
notb = find(notb);
Q11 = Q(notb, notb);
Q12 = Q(notb, b);

for j=1:length(b)
    % #H x 1 = (a x a) \ (axb * bx1) = a x 1
    W(notb,j) = Q11 \ -(Q12*bc(:,j));
end
    
end