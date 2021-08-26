clear all
clc
close all

model = createpde('structural','static-planestress');

%% Use our own mesh
meshPath = '../data/woody.obj';
[V,F] = readOBJ(meshPath);

% figure;
% fig = tsurf(F,V);
% axis equal;

V = V(:,1:2); % the input mesh contains redundant zero third column
geometryFromMesh(model,V',F');

figure;
% pdeplot(model)
pdegplot(model,'EdgeLabels','on','VertexLabels','on');
ylim([-1,11])
axis equal
title 'Geometry With Edge Labels Displayed';
hold on

%% Model parameters
structuralProperties(model,'YoungsModulus',200E3,'PoissonsRatio',0.25);

%% Loads
structuralBoundaryLoad(model,'Edge',6,'SurfaceTraction',[0;-10]);

%% Supports
% structuralBC(model,'Edge',1:9,'XDisplacement',0);
% structuralBC(model,'Edge',1:9,'YDisplacement',0);
structuralBC(model,'Edge',[8],'Constraint','fixed');
% structuralBC(model,'Edge',[8],'Constraint','fixed');

%% Plot the FEM mesh
figure
pdemesh(model);

R = solve(model);

%% Plot results
% figure
% pdeplot(model,'XYData',R.Stress.sxx,'ColorMap','jet')
% axis equal
% title 'Normal Stress Along x-Direction';
figure
pdeplot(model,'XYData',R.Displacement.Magnitude,'ColorMap','jet')
axis equal
title 'Deformation magnitude';

figure
pdeplot(model,'XYData',R.Displacement.ux,'ColorMap','jet')
axis equal
title 'Deformation ux';

figure
pdeplot(model,'XYData',R.Displacement.uy,'ColorMap','jet')
axis equal
title 'Deformation uy';