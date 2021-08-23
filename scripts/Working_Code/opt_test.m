% Finds optimal shape deformed using selected handles based on complexity, i.e. finds the most complex deformation

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

% geet handle indices b
[~,b] = min(D);

% TODO: (unbounded) biharmonic weights
W = compute_skinning_weight(V,F,b);


%A = [eye(length(C)*3); eye(length(C)*3)*-1];
%b = [ones(length(C)*2,1)*25; ones(length(C),1)*(pi/4); ones(length(C)*2,1)*25; ones(length(C),1)*(pi/4)];
lb = [ones(length(C)*2,1)*-100; ones(length(C),1)*(-pi/4)];
ub = [ones(length(C)*2,1)*100; ones(length(C),1)*(pi/4)];
x0 = zeros(length(C)*3,1);
%[xsol, fval] = fmincon(@f, x0, [], [], [], [], lb, ub);
[xsol,fval,history,searchdir] = runfmincon(x0, lb, ub);


new_C = C + [xsol(1:length(C)), xsol((length(C)+1):(length(C)*2))];
TR = skinning_transformations(C, P, [], new_C, xsol((2*length(C)+1):end));
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
TR = skinning_transformations(C, P, [], new_C, X((2*length(C)+1):end));
new_V = linear_blend_skinning(V(:,1:2), TR, W);

c = complexity(new_V, F);

end

function [c, ceq] = nonlcon(X)

global C;
global P;
global V;
global F;
global W;

new_C = C + [X(1:length(C)), X((length(C)+1):(length(C)*2))];
TR = skinning_transformations(C, P, [], new_C, X((2*length(C)+1):end));
new_V = linear_blend_skinning(V(:,1:2), TR, W);

c = 0.95 - acceptability(new_V, F);
ceq = [];

end

function [xsol,fval,history,searchdir] = runfmincon(x0, lb, ub)
 
% Set up shared variables with outfun
history.x = [];
history.fval = [];
searchdir = [];
 
% Call optimization
options = optimoptions(@fmincon,'OutputFcn',@outfun);
[xsol,fval] = fmincon(@f,x0,[],[],[],[],lb,ub,@nonlcon,options);
 
 function stop = outfun(x,optimValues,state)
     global V;
     global C;
     global P;
     global W;
     global F;
     
     stop = false;
 
     switch state
         case 'init'
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
            history.fval = [history.fval, optimValues.fval];
            history.x = [history.x, x];
         % Plot deformed mesh
            new_C = C + [x(1:length(C)), x((length(C)+1):(length(C)*2))];
            TR = skinning_transformations(C, P, [], new_C, x((2*length(C)+1):end));
            new_V = linear_blend_skinning(V(:,1:2), TR, W);
            tsurf(F, new_V); axis equal;
            hold on;
            scatter3( new_C(:,1),new_C(:,2),0.1+0*new_C(:,1), 'o',...
                'MarkerFaceColor', [0.9 0.8 0.1], 'MarkerEdgeColor','k',...
                'LineWidth',2,'SizeData',100);
            drawnow;
            hold off;
         case 'done'
         otherwise
     end
 end

end
