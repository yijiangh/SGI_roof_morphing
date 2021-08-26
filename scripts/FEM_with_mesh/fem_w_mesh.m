
clear all
clc
close all

numberOfPDE = 2;
model = createpde(numberOfPDE);

%% Problem Parameters
%%
E = 1.0e6; % modulus of elasticity
nu = .3; % Poisson's ratio
thick = .1; % plate thickness
len = 10.0; % side length for the square plate
hmax = len/20; % mesh size parameter
D = E*thick^3/(12*(1 - nu^2));
pres = 2; % external pressure

%% Geometry Creation
% For a single square, the geometry and mesh are easily defined as shown below.

% gdm = [3 4 0 len len 0 0 0 len len]';
% g = decsg(gdm,'S1',('S1')');
% Create a geometry entity.
% geometryFromEdges(model,g);

%% Use our own mesh
meshPath = '../data/woody.obj';
[V,F] = readOBJ(meshPath);

% figure;
% fig = tsurf(F,V);
% axis equal;

V = V(:,1:2); % the input mesh contains redundant zero third column
geometryFromMesh(model,V',F');

%% 
% Plot the geometry and display the edge labels for use in the boundary 
% condition definition.

figure;
% pdeplot(model)
pdegplot(model,'EdgeLabels','on','VertexLabels','on');
ylim([-1,11])
axis equal
title 'Geometry With Edge Labels Displayed';

c = [1 0 1 D 0 D]';
a = [0 0 1 0]';
f = [0 pres]';
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',f);

%% Boundary Conditions
%%
k = 1e7; % spring stiffness
%% 
% Define distributed springs on all four edges.

% 1:10 is the edge label shown in the figure above
% bOuter = applyBoundaryCondition(model,'neumann','Edge',[2,3,9],...
%                                      'g',[0 0],'q',[0 0; k 0]);

bOuter = applyBoundaryCondition(model,'neumann','Vertex',1:2,...
                                     'g',[0 0],'q',[0 0; k 0]);
                                 
% bOuter = applyBoundaryCondition(model,'dirichlet','Edge',(1:4),...
%                                      'u',[0 0]);

%% Mesh generation
%%
% we bypass the meshing part
% tic
% generateMesh(model, 'Hmax', hmax);
% toc

%%

tic
res = solvepde(model);
toc

u = res.NodalSolution;
numNodes = size(model.Mesh.Nodes,2);
figure
pdeplot(model,'XYData',u(1:numNodes),'Contour','on');
title 'Transverse Deflection'

numNodes = size(model.Mesh.Nodes,2);
fprintf('Transverse deflection at plate center(PDE Toolbox) = %12.4e\n',...
                                                  min(u(1:numNodes,1)));
%% 
% Compute analytical solution.

% wMax = -.0138*pres*len^4/(E*thick^3);
% fprintf('Transverse deflection at plate center(analytical) = %12.4e\n', wMax);
%% 
% _Copyright 2012-2015 The MathWorks, Inc._