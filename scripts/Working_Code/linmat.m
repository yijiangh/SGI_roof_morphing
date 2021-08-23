function [A, b] = linmat(V, W)
%LINMAT returns A and b from the linear equation A x + b = V_new which computes the locations of the new vertices after the handles are translated by x. Both x and V_new are one dimensional vectors with x coordinates stored in the first half of the vector and y coordinates stored in the second half

[m, n] = size(W);

A = zeros(m, 3*n);
A(:,1:3:end) = W .* V(:,1);
A(:,2:3:end) = W .* V(:,2);
A(:,3:3:end) = W;

bigA = [A, zeros(m,3*n); zeros(m,3*n), A];

Ta = zeros(6*n,2*n);
for i=1:2*n
    Ta(3*i,i) = 1;
end

Tb = zeros(6*n,1);
for i=1:n
    Tb(3*i-2) = 1;
    Tb(3*n + 3*i-1) = 1;
end

A = bigA * Ta;
b = bigA * Tb;

end

