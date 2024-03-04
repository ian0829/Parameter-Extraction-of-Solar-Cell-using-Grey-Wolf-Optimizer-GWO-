
% PSO practice

%% ---A--- Set up scenario
x1 = -3:0.1:3;
x2 = -3:0.1:3;
[x1,x2] = meshgrid(x1, x2);

y = (4-2.1.*x1.^2+(x1.^4)/3) .* x1.^2+x1.*x2+ (-4+4.*x2.^2) .* x2.^2;

figure;
surf(x1, x2, y);
hold on;
scatter3 (0.0898, -0.7126, -1.0316, 120, 'MarkerFaceColor', 'y');
scatter3 (-0.0898, 0.7126, -1.0316, 120, 'MarkerFaceColor', 'y');

%% --- B--- PSO function applies
%【1】 formulate the problem
f = @(x1,x2) (4-2.1.*x1.^2+(x1.^4)/3) .* x1.^2+x1.*x2+ (-4+4.*x2.^2) .* x2.^2;
fun = @(x)f(x(1),x(2));

%【2】 define lower bound 
lb = [-3; -3];
%【3】 define upper bound 
ub = -lb;
%【4】 number of variables
nvars = 2;
%【5】  use function "particleswarm"
[x, fval, exitflag] = particleswarm (fun, nvars, lb, ub);
%【6】  output the best result as minFitness
minFitness = fval;
%【7】 output the position of the best result as minPosition
minPosition = x;

% mark the best point
scatter3(minPosition(1),minPosition(2), minFitness ,25,'MarkerEdgeColor','r','MarkerFaceColor' ,'r');
hold on;
view(-113,60);

