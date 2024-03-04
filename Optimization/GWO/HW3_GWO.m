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
scatter3 (0.0898, -0.7126, -1.0316, 120, 'MarkerFaceColor', 'y');
scatter3 (-0.0898, 0.7126, -1.0316, 120, 'MarkerFaceColor', 'y');
hold on;
view(82, 53);

%% --------------------------------------
% ----------- GWO Parameters ------------
%----------------------------------------
NE              = 2;   % number of elements
N               = 20; % number of particles

iter            = 100; % number of iterations
iterCount       = 0;   % counter of iterations
a_max           = 2; 
a_min           = 0; 

%  -------------------------------------------
%  ------------ Boundary Settings ------------
%  -------------------------------------------
P_Max           = [3, 3];       % [ upper bound of the searching space ]
P_Min           = [-3, -3];     % [ lower bound of the searching space ]
P               = zeros (N, NE); 
Alpha_pos       = zeros (iter, NE);  %Location
Beta_pos        = zeros (iter, NE);  %Location
Delta_pos       = zeros (iter, NE);  %Location
Alpha_score     = zeros (iter, 1);
Beta_score      = zeros (iter, 1);
Delta_score     = zeros (iter, 1);
Fitness         = zeros (N, 1);      %Score

%-----------------------------------------
% ------------ Initialization-------------
%-----------------------------------------
 for i = 1 : N
     for j = 1 : NE
        P(i, j) = P_Min(1,j) + (P_Max(1,j) - P_Min(1,j)) * rand();
     end
 end
     
 while ( iterCount < iter )
%-------------------------------------
% ------------ Evaluation ------------
%-------------------------------------
for i = 1 : N
    Fitness(i, 1) = camel6(P(i,:));
end
move = scatter3(P(:, 1), P(:,2), Fitness, 60, 'MarkerFaceColor', 'r');
pause(0.5);
delete(move);

% ----------------------------------------------------------------------
%               Update each search agent location & Fitness value
% ----------------------------------------------------------------------
Seq = Fitness; 
for i = 1 : N                   % Bubble sort
    for j = 1 : N-1
        if Seq(j,1) > Seq(j+1,1)
            temp = Seq(j,1);
            Seq(j) = Seq(j+1,1);
            Seq(j+1,1) = temp;
        end
    end
end
for i = 1 : N
    if Fitness(i,1) == Seq (1,1)
        Alpha_score(iterCount+1, 1) = Fitness(i, 1); % Update alpha
        Alpha_pos(iterCount+1, :) = P(i,:);
    elseif Fitness(i,1) == Seq (2,1)
        Beta_score(iterCount+1, 1) = Fitness(i,1); % Update beta
        Beta_pos(iterCount+1, :) = P(i,:);
    elseif Fitness(i,1) == Seq (3,1)
        Delta_score(iterCount+1, 1) = Fitness(i, 1); % Update delta
        Delta_pos(iterCount+1, :) = P(i,:);
    end
end

% ---------------------------------------------------
%               Update locations
% ---------------------------------------------------
for i = 1 : N
    a = a_max - (iterCount/iter) * (a_max-a_min);

    A1 = 2.*a.*rand(1,2)-a;
    C1 = 2.*rand(1,2);
    D_alpha = abs(C1.*Alpha_pos(iterCount+1,:)- P(i,:));
    X1 = Alpha_pos(iterCount+1,:) - A1.*D_alpha;
    
    A2 = 2.*a.*rand(1,2)-a;
    C2 = 2.*rand(1,2);    
    D_beta = abs(C2.*Beta_pos(iterCount+1,:)- P(i,:));
    X2 = Beta_pos(iterCount+1,:) - A2.*D_beta;
    
    A3 = 2.*a.*rand(1,2)-a;
    C3 = 2.*rand(1,2);    
    D_delta = abs(C3.*Delta_pos(iterCount+1,:)- P(i,:));
    X3 = Delta_pos(iterCount+1,:) - A3.*D_delta;
   
    P(i, :) = (X1+X2+X3)/3;
    
    P(i, :) = max(P(i, :), P_Min(1,:));
    P(i, :) = min(P(i, :), P_Max(1,:));
end

% ---------------------------------------------------
%               Next Iteration
% ---------------------------------------------------
% display iteration (iterCount+1)

iterCount = iterCount + 1;    

    
 end
scatter3(P(:, 1), P(:,2), Fitness, 60, 'MarkerFaceColor', 'r');
hold;
%% Plot the convergence history
figure;
plot(Alpha_score, 'b');
grid on;
hold on;
xlabel('Generation');
ylabel('Fitness');
title('Convergence History');

toc  