clc; clear all; close all;

%meshPath = '../data/woody.obj';
%[V,F] = readOBJ(meshPath);
%V = V(:,1:2); % the input mesh contains redundant zero third column

pts = get_pencil_curves(0.0005);
E = [(1:size(pts,1))',[(2:size(pts,1))';1]];
[V,F] = triangle(pts,E, [], 'MaxArea', 0.0005, 'Flags','-q20');

% get handles
fig = tsurf(F,V);
axis equal;
fprintf( ...
    ['Point Handle Selection: \n' ...
    '- CLICK the mesh to add point handls \n', ...
    '- BACKSPACE to remvoe the previous selection\n', ... 
    '- ENTER to finish selection\n'] ...
    );
try
  [Cx,Cy] = getpts;
catch e
  return  % quit early, stop script
end

C = [Cx,Cy]; % handle locations 

% compute pairwise distance
D = zeros(size(V,1), size(C,1));
for ii = 1:size(C,1)
    D(:,ii) = sqrt(sum((V - C(ii,:)).^2,2));
end

% get handle indices b
[~,b] = min(D);

% TODO: (unbounded) biharmonic weights
W = compute_skinning_weight(V,F,b);

% TODO: linear blended skinning
DeformOBJ = deform_GUI(V,F,C,W);
