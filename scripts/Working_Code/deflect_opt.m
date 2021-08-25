% Finds optimal shape deformed using selected handles based on deflection
% while constraining AREA

clc; clear; close all;

global C;
global P;
global V;
global F;
global W;

% Get V and F
[V, F] = readOBJ('woody.obj');
V = V(:,1:2);
V = (V - mean(V)) * 0.2;

% Get handles (automate this in future)
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
W = compute_skinning_weight(V,F,b);

% Set up optimization constraints and run fmincon
lb = ones(length(C)*2,1)*-10;
ub = ones(length(C)*2,1)*10;
x0 = zeros(length(C)*2,1);
[xsol,fval,history,searchdir] = runfmincon(x0, lb, ub);

% Compute new mesh and values after optimization is complete
new_C = C + [xsol(1:length(C)), xsol((length(C)+1):end)];
TR = skinning_transformations(C, P, [], new_C, zeros(length(C),1));
new_V = linear_blend_skinning(V(:,1:2), TR, W);

[model, u] = deflection(new_V, F);


% Optimization function
function d = f(X)

global C;
global P;
global V;
global F;
global W;

% Compute new mesh vertices
new_C = C + [X(1:length(C)), X((length(C)+1):end)];
TR = skinning_transformations(C, P, [], new_C, zeros(length(C),1));
new_V = linear_blend_skinning(V(:,1:2), TR, W);

try
    % Compute and return max deflection
    [model, u] = deflection(new_V, F);
    numNodes = size(model.Mesh.Nodes,2);
    d = -min(u(1:numNodes,1));
catch
    d = 99999; % Catch if triangles are too smooshed. Return large number
end
end

% Nonlinear area constraint
function [c, ceq] = nonlcon(X)

global C;
global P;
global V;
global F;
global W;

new_C = C + [X(1:length(C)), X((length(C)+1):end)];
TR = skinning_transformations(C, P, [], new_C, zeros(length(C),1));
new_V = linear_blend_skinning(V(:,1:2), TR, W);

c = (sum(doublearea(new_V, F)) / 2) - (sum(doublearea(V, F)) / 2) + 25;
ceq = [];

end

function [xsol,fval,history,searchdir] = runfmincon(x0, lb, ub)
 
% Set up shared variables with outfun
history.x = [];
history.fval = [];
searchdir = [];
 
% Call optimization
options = optimoptions(@fmincon,'OutputFcn',@outfun);
[xsol,fval] = fmincon(@f,x0,[],[],[],[],lb,ub,@nonlcon, options);
 
 function stop = outfun(x,optimValues,state)
     global V;
     global C;
     global P;
     global W;
     global F;
     global i;
     
     stop = false;
 
     switch state
         case 'init'
             i = 0;
         case 'iter'
             i = i+1;
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
            history.fval = [history.fval, optimValues.fval];
            history.x = [history.x, x];
         % Plot deformed mesh
            new_C = C + [x(1:length(C)), x((length(C)+1):end)];
            TR = skinning_transformations(C,P,[],new_C,zeros(length(C),1));
            new_V = linear_blend_skinning(V(:,1:2), TR, W);
            subplot(131);
            tsurf(F, new_V); axis equal;
            hold on;
            scatter3( new_C(:,1),new_C(:,2),0.1+0*new_C(:,1), 'o',...
                'MarkerFaceColor', [0.9 0.8 0.1], 'MarkerEdgeColor','k',...
                'LineWidth',2,'SizeData',100);
            drawnow;
            hold off;
            
            [model, u] = deflection(new_V, F);
            subplot(132);
            pdeplot(model, 'XYdata', u, 'Zdata', u);
            
            subplot(133);
            plot((1:i), history.fval, 'b.')
         case 'done'
         otherwise
     end
 end

end
