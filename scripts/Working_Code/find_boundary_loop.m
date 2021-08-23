function [boundary_loop] = find_boundary_loop(V, F)
%give ordered list of points along boundary of mesh
O = outline(F);

boundary_loop = zeros(length(O), 1);
v = O(1,1);
for i=1:length(boundary_loop)
    boundary_loop(i) = v;
    ind = find(O(:,1) == v);
    v = O(ind,2);
end

boundary_loop = V(boundary_loop, :);

end