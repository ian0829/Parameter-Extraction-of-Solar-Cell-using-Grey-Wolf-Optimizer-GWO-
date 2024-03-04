clear;
clc;
tic;
%-------------------------------------------
% ----------- Import Measured Data----------
%--------------------------------------------
load('IPV.mat');
load('VPV.mat');
global Vpv0 Ipv0
Vpv0 = VPV; % measured voltage data 
Ipv0 = IPV; % measured current data


plot(Vpv0,Ipv0, 'r', 'LineWidth', 3) 
legend('measured')
xlim([0 70])
ylim([0 10])
ylabel('I [A]')
xlabel('U [V]')
title('Solar Cell I-V Characteristic')
hold on;
%% --------------------------------------------------
% ----------- GWO Parameters ------------
%----------------------------------------
NE              = 5;   % number of elements
N               = 100;  % number of wolfs

iter            = 1000; % number of iterations
iterCount       = 0;   % counter of iterations
a_max           = 2; 
a_min           = 0; 

% ------------------------------------------
% ----------- Boundary Settings ------------
%-------------------------------------------
%P = [Rs, Rsh, Iph, Is, n];
P_Max = [  1, 300, 10, 0.0001, 2];
P_Min = [0.1, 288,  5,      0, 1]; 
P               = zeros (N, NE);
Alpha_pos       = zeros (iter, NE);  %Location
Beta_pos        = zeros (iter, NE);  %Location
Delta_pos       = zeros (iter, NE);  %Location
Alpha_score     = zeros (iter, 1);
Beta_score      = zeros (iter, 1);
Delta_score     = zeros (iter, 1);
rmse            = zeros (N, 1);

%-----------------------------------------
% ------------ Initialization-------------
%-----------------------------------------
 for i = 1 : N
     for j = 1 : NE
        P(i, j) = (P_Min(1,j) + (P_Max(1,j) - P_Min(1,j)) * round(rand(),9));
     end
 end
 while ( iterCount < iter )
%-------------------------------------
% ------------ Evaluation ------------
%-------------------------------------
for i = 1 : N
    [rmse(i,1), Ipv(i,:)] = fobj(P(i,:));
end
%{
move = plot(Vpv0,Ipv, 'g');
legend('measured','data')
pause(0.1);
delete(move);
%}

% ----------------------------------------------------------------------
%               Update each search agent location & Fitness value
% ----------------------------------------------------------------------
Seq = sort(rmse, 'ascend');
%{
% Bubble sort
for i = 1 : N
    for j = 1 : N-1
        if Seq(j,1) > Seq(j+1,1)
            temp = Seq(j,1);
            Seq(j) = Seq(j+1,1);
            Seq(j+1,1) = temp;
        end
    end
end
%}
for i = 1 : N
    if rmse(i,1) == Seq (1,1)
        Alpha_score(iterCount+1, 1) = rmse(i, 1); % Update alpha
        Alpha_pos(iterCount+1, :) = P(i,:);
    elseif rmse(i,1) == Seq (2,1)
        Beta_score(iterCount+1, 1) = rmse(i,1); % Update beta
        Beta_pos(iterCount+1, :) = P(i,:);
    elseif rmse(i,1) == Seq (3,1)
        Delta_score(iterCount+1, 1) = rmse(i, 1); % Update delta
        Delta_pos(iterCount+1, :) = P(i,:);
    end
end

% ---------------------------------------------------
%               Update locations
% ---------------------------------------------------
for i = 1 : N
    a = a_max - (iterCount/iter) * (a_max-a_min);

    A1 = 2.*a.*round(rand(1,NE), 9)-a;
    C1 = 2.*round(rand(1,NE), 9);
    D_alpha = abs(C1.*Alpha_pos(iterCount+1,:)- P(i,:));
    X1 = Alpha_pos(iterCount+1,:) - A1.*D_alpha;
    
    A2 = 2.*a.*round(rand(1,NE), 9)-a;
    C2 = 2.*round(rand(1,NE), 9);    
    D_beta = abs(C2.*Beta_pos(iterCount+1,:)- P(i,:));
    X2 = Beta_pos(iterCount+1,:) - A2.*D_beta;
    
    A3 = 2.*a.*round(rand(1,NE), 9)-a;
    C3 = 2.*round(rand(1,NE), 9);    
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
plot(Vpv0,Ipv, 'b');
legend('measured','data');
hold;
%% Plot the convergence history
figure;
plot(Alpha_score, 'b');
grid on;
hold on;
xlabel('Generation');
ylabel('RMSE');
title('Convergence History');

toc;  