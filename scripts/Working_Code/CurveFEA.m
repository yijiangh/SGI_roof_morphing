
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MIT Summer Geometry Institute 2021
%   Project 003: Design Optimization Using Shape Morphing
%   Mentors: Caitlin Mueller & Yijiang
%
%   Finite Element Analysis [User-Drawn Curve]
%   
%   Author(s): Marcus Vidaurri, Lily Kimble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Clear Cache

clear all %#ok<*CLALL>
close all
clc
rng('shuffle')

fprintf('~ Finite Element Analysis for a User-Drawn Curve ~\n\n')

%====================================================================
%%  FEA: User-Drawn Curve

% NOTE: Make sure to draw an open curve (the FEA will close it for meshing)

% Draw Curve
[V,Edges,cid] = get_pencil_curves(1e-6); % Type 'n' and Enter to finish


% Create Polygon
id = 2;                         % shape ID: 2 = polygon of N sides
N = size(V,1);                  % number of sides from boundary edges
x = V(:,1)';                    % x-coords of boundary vertices
y = V(:,2)';                    % y-coords of boundary vertices
gdm = [id, N, x, y]';           % geometry description matrix
dgm = decsg(gdm,'S1',('S1')');  % decomposed geometry matrix

model = createpde(2);           % setup PDE
geometryFromEdges(model,dgm);   % add the geometry representation to PDE

figure(2)
pdegplot(model,'EdgeLabels','off') % Change 'off' to 'on' to show labels
% xlim([min(x)-1 max(x)+1])
% ylim([min(y)-1 max(y)+1])
axis equal
title 'Geometry With Edge Labels Displayed'

E = 1.0e6; % Modulus of elasticity
nu = 0.3; % Poisson's ratio
thick = 0.1; % Plate thickness
pres = 500; % External pressure

D = E*thick^3/(12*(1 - nu^2));

c = [1 0 1 D 0 D]';
a = [0 0 1 0]';
f = [0 pres]';
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',f);

k = 1e7; % spring stiffness

% Boundary Conditions
bOuter = applyBoundaryCondition(model,'neumann','Edge',(1:N),...
                                                 'g',[0 0],'q',[0 0; k 0]);

generateMesh(model);

res = solvepde(model);

u = res.NodalSolution;

numNodes = size(model.Mesh.Nodes,2);

figure(3)

subplot(2,1,1)
pdeplot(model,'XYData',u(:,1),'Contour','on') 
title 'Transverse Deflection'

subplot(2,1,2)
pdeplot(model,'XYData',u,'Zdata',u)
title '3D deflection'

wMax_n = min(u(1:numNodes,1));  % max deflection [numerical solution]

fprintf(['\nMax Deflection: ' num2str(wMax_n)])


%====================================================================
%%  End of File

fprintf('\n\nEnd of File. All plots generated.\n\n')
