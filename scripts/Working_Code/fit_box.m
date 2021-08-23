% Deforms shape using selected handles to fill as much of a box without
% exceeding box boundaries. Attempts to avoid overlap and extreem
% distortion by monitering aspect ratio of triangles

clc; clear all; close all;

global C;
global P;
global V;
global F;
global W;


[V, F] = readOBJ('woody.obj');


V = V(:,1:2);
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
P = 1:size(C,1);

% compute pairwise distance
D = zeros(size(V,1), size(C,1));
for ii = 1:size(C,1)
    D(:,ii) = sqrt(sum((V - C(ii,:)).^2,2));
end

% get handle indices b
[~,b] = min(D);

% TODO: (unbounded) biharmonic weights
W = compute_skinning_weight(V,F,b);


x0 = zeros(length(C)*2,1);
[A, b] = linmat(V, W);
A = [A;-1*A];
b = [-1*b+410; b + 10];
%[xsol, fval] = fmincon(@f, x0, [], [], [], [], lb, ub);
[xsol,fval,history,searchdir] = runfmincon(x0, A, b);


new_C = C + [xsol(1:length(C)), xsol((length(C)+1):(length(C)*2))];
TR = skinning_transformations(C, P, [], new_C, zeros(length(C),1));
new_V = linear_blend_skinning(V(:,1:2), TR, W);

c = complexity(new_V, F);

figure; subplot(121);
tsurf(F, V); axis equal;
hold  on;
scatter3( C(:,1),C(:,2),0.1+0*C(:,1), 'o',...
    'MarkerFaceColor', [0.9 0.8 0.1], 'MarkerEdgeColor','k',...
    'LineWidth',2,'SizeData',100);
hold off;

subplot(122);
tsurf(F, new_V); axis equal;
hold on;
scatter3( new_C(:,1),new_C(:,2),0.1+0*new_C(:,1), 'o',...
    'MarkerFaceColor', [0.9 0.8 0.1], 'MarkerEdgeColor','k',...
    'LineWidth',2,'SizeData',100);
hold off;

function c = f(X)

global C;
global P;
global V;
global F;
global W;

new_C = C + [X(1:length(C)), X((length(C)+1):(length(C)*2))];
TR = skinning_transformations(C, P, [], new_C, zeros(length(C),1));
new_V = linear_blend_skinning(V(:,1:2), TR, W);

c = -sum(doublearea(new_V, F)) / 2;

c = c * (acceptability(new_V, F))^2;

end

function [c, ceq] = nonlcon(X)

global C;
global P;
global V;
global F;
global W;

new_C = C + [X(1:length(C)), X((length(C)+1):(length(C)*2))];
TR = skinning_transformations(C, P, [], new_C, zeros(length(C),1));
new_V = linear_blend_skinning(V(:,1:2), TR, W);

c = 0.95 - acceptability(new_V, F);
ceq = [];

end

function [xsol,fval,history,searchdir] = runfmincon(x0, A, b)
 
% Set up shared variables with outfun
history.x = [];
history.fval = [];
searchdir = [];
 
% Call optimization
options = optimoptions(@fmincon,'OutputFcn',@outfun);
[xsol,fval] = fmincon(@f,x0,A,b,[],[],[],[],[],options);
 
 function stop = outfun(x,optimValues,state)
     global V;
     global C;
     global P;
     global W;
     global F;
     
     stop = false;
 
     switch state
         case 'init'
             %figure;
             %hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval, optimValues.fval];
           history.x = [history.x, x];
         % Concatenate current search direction with 
         % searchdir.
           %searchdir = [searchdir;... 
           %             optimValues.searchdirection'];
           %plot(x(1),x(2),'o');
         % Label points with iteration number and add title.
         % Add .15 to x(1) to separate label from plotted 'o'.
           %text(x(1)+.15,x(2),... 
           %     num2str(optimValues.iteration));
           %title('Sequence of Points Computed by fmincon');
            new_C = C + [x(1:length(C)), x((length(C)+1):(length(C)*2))];
            TR = skinning_transformations(C, P, [], new_C, zeros(length(C),1));
            new_V = linear_blend_skinning(V(:,1:2), TR, W);
            tsurf(F, new_V); axis equal;
            hold on;
            scatter3( new_C(:,1),new_C(:,2),0.1+0*new_C(:,1), 'o',...
                'MarkerFaceColor', [0.9 0.8 0.1], 'MarkerEdgeColor','k',...
                'LineWidth',2,'SizeData',100);
            drawnow;
            hold off;
         case 'done'
             %hold off
         otherwise
     end
 end

end