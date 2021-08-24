function [model, u] = deflection(V, F)
%DEFLECTION Summary of this function goes here
%   Detailed explanation goes here

model = createpde(2);

pts = find_boundary_loop(V, F);

gdm = [2, length(pts), pts(:,1)', pts(:,2)']';

g = decsg(gdm,'S1', ('S1')');
geometryFromEdges(model,g);

E = 1.0e6; % Modulus of elasticity
nu = 0.3; % Poisson's ratio
thick = 0.1; % Plate thickness
pres = 2; % External pressure

D = E*thick^3/(12*(1 - nu^2));

c = [1 0 1 D 0 D]';
a = [0 0 1 0]';
f = [0 pres]';
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',f);

k = 1e7;

bOuter = applyBoundaryCondition(model,'neumann','Edge',(1:length(pts)),...
                                     'g',[0 0],'q',[0 0; k 0]);
                                 
generateMesh(model);

res = solvepde(model);

u = res.NodalSolution;

numNodes = size(model.Mesh.Nodes,2);
%pdeplot(model,'XYData',u(:,1),'Contour','on')
%drawnow();

numNodes = size(model.Mesh.Nodes,2);
wMax = min(u(1:numNodes,1));

end

