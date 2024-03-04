clear
clc
tic
% typical PSO
% evaluated with test function

%% ---------------------------------------------------
% ----------- Scenario setup ------------
%----------------------------------------
x1 = -3:0.1:3;
x2 = -3:0.1:3;
[x1,x2] = meshgrid(x1, x2);

y = (4-2.1.*x1.^2+(x1.^4)/3) .* x1.^2+x1.*x2+ (-4+4.*x2.^2) .* x2.^2;

figure;
surf(x1, x2, y);
hold on;
scatter3 (0.0898, -0.7126, -1.0316, 120, 'MarkerFaceColor', 'y'); %用黃點繪出最佳解座標
scatter3 (-0.0898, 0.7126, -1.0316, 120, 'MarkerFaceColor', 'y'); %用黃點繪出最佳解座標
hold on;
view(82, 53);

%% --------------------------------------
% ----------- PSO Parameters ------------
%----------------------------------------
NE              = 2;   % number of elements
N               = 100; % number of particles

iter            = 100; % number of iterations
iterCount       = 0;   % counter of iterations
c1              = 1; 
c2              = 2; 
w_max           = 1; 
w_min           = 0.1; 

%  -------------------------------------------
%  ------------ Boundary Settings ------------
%  -------------------------------------------
P_Max           = [3, 3];       % [ upper bound of the searching space ]
P_Min           = [-3, -3];     % [ lower bound of the searching space ]
V_Max           = (P_Max - P_Min) / 10;  % [ max movement ]
V_Min           = -V_Max;                % [ min movement ]

P               = zeros (N, NE); 
V               = zeros (N, NE);
pBest           = zeros (N, NE);     %Location
gBest           = zeros (iter, NE);  %Location
Fitness         = zeros (N, 1);      %Score
Fitness_pBest   = zeros (N, 1);      %Score
Fitness_gBest   = zeros (iter, 1);   %Score

%-----------------------------------------
% ------------ Initialization-------------
%-----------------------------------------
 for i = 1 : N
     for j = 1 : NE
        P(i, j) = P_Min(1,j) + (P_Max(1,j) - P_Min(1,j)) * rand();
        V(i, j) = V_Min(1,j) + (V_Max(1,j) - V_Min(1,j)) * rand();
     end
 end
     
 while ( iterCount < iter )
%-------------------------------------
% ------------ Evaluation ------------
%-------------------------------------
for i = 1 : N
    Fitness (i, 1) = camel6(P(i,:));
end
move = scatter3(P(:, 1), P(:,2), Fitness, 60, 'MarkerFaceColor', 'r'); %將每一次迭代的粒子位置用紅點繪出
pause(0.5); %顯示0.5秒
delete(move); %刪除當次迭代標記

% update particle location



% ----------------------------------------------------------------------
%               Update pBest and pBest location & Fitness value
% ----------------------------------------------------------------------
for i = 1 : N
    if Fitness_pBest(i, 1) > Fitness(i, 1)
       Fitness_pBest(i, 1) = Fitness(i, 1);
       pBest(i, :) = P(i, :); 
    end
    if Fitness_gBest(iterCount+1, 1) > Fitness_pBest (i, 1)
       Fitness_gBest(iterCount+1, 1) = Fitness_pBest (i, 1);
       gBest(iterCount+1, :) = pBest(i, :);
    end
end

% ---------------------------------------------------------
%              Calculation velocities 
% ---------------------------------------------------------%
w = w_max - (iterCount/iter) * (w_max-w_min);
for i = 1 : N
    V(i, :) = w*V(i, :) + rand()*c1.*(pBest(i, :)-P(i, :)) + rand()*c2.*(gBest(iterCount+1, :)-P(i, :));
    V(i, :) = max(V(i, :), V_Min(1,:));
    V(i, :) = min(V(i, :), V_Max(1,:));
end

% ---------------------------------------------------
%               Update locations
% ---------------------------------------------------
for i = 1 : N
    P(i, :) = P(i, :) + V(i, :);
    P(i, :) = max(P(i, :), P_Min(1,:));
    P(i, :) = min(P(i, :), P_Max(1,:));
end



% ---------------------------------------------------
%               Next Iteration
% ---------------------------------------------------
% display iteration (iterCount+1)
% display best value (Fitness_gBest)
% display best setting (gBest)


iterCount = iterCount + 1;    

    
 end
scatter3(P(:, 1), P(:,2), Fitness, 60, 'MarkerFaceColor', 'r');
hold;
%% Plot the convergence history
figure;
plot(Fitness_gBest, 'b');   %標記群體最佳的分數
grid on;
hold on;
xlabel('Generation');
ylabel('Fitness');
title('Convergence History');

toc  