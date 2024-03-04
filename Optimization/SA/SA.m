% Simulated Anealing algorithm
% 2017/11/30
clear;
clc;
tic;
%%
figure; 
u = -10:0.1:10;
for i  =  1:size(u,2)
    y(i) = levy(u(i));
end

plot(u,y);
hold on;
title('Simulated Anealing (SA) for levy function');
xlabel('x');
ylabel('y');

%% --------------------------------------
% ----------- SA Parameters ------------
%----------------------------------------
%T           = 1000;
NE          = 1;                % number of element
x_Max       = 10*ones(1,NE);
x_Min       = -10*ones(1,NE);
v_Max       = 5*ones(1,NE);
v_Min       = -5*ones(1,NE);
Fitness     = 0;
iterCount   = 0;
iter        = 500;
NEWx        = zeros( iter , NE);
t_init      = 500;


while ( iterCount < iter )

    if  iterCount == 0
        %-----------------------------------------
        % ------------ Initialization-------------
        %-----------------------------------------   
        for i = 1 : NE
            x( iterCount+1, i) = x_Min(i) + (x_Max(i)-x_Min(i))*rand();
        end
        Fitness( iterCount+1 ) = levy( x( iterCount+1, :) );
         T = t_init*(0.9^iterCount);
        move = scatter(x(:,1), Fitness( iterCount+1 ) ,30,'*','r','LineWidth',2.5);
        pause(0.1);
        delete(move);


    else

        %-------------------------------------
        % ----------- Evaluation ------------
        %--------------------------------------

        FitCandidate( iterCount ) = levy( NEWx(iterCount,:) );

        move = scatter(NEWx(iterCount,:) , FitCandidate( iterCount ) ,30,'*','r','LineWidth',2.5);
        pause(0.1);
        delete(move);


        deltaE = FitCandidate( iterCount ) - Fitness( iterCount );

        if deltaE < 0
            Fitness( iterCount+1 ) = FitCandidate( iterCount );
            x(iterCount+1,:) = NEWx(iterCount,:);
        else
            if rand() < exp(-deltaE/T)
                Fitness( iterCount+1 ) = FitCandidate( iterCount );
                x(iterCount+1,:) = NEWx(iterCount,:);
            else
                Fitness( iterCount+1 ) = Fitness( iterCount );
                x(iterCount+1,:) = x(iterCount,:);
            end
            T = t_init*(0.9^iterCount);
        end
    end


    %----------------------------------------------
    % ------------  Movement   Update -------------
    %----------------------------------------------

    for i = 1 : NE
        v( iterCount+1 , i) = v_Min(i) + (v_Max(i)-v_Min(i))*rand();
        NEWx(iterCount+1,i) = x(iterCount+1,i) + v(iterCount+1,i);
    end
    
    NEWx(iterCount+1,:) = min(x_Max, NEWx(iterCount+1,i));
    NEWx(iterCount+1,:) = max(x_Min, NEWx(iterCount+1,i));
    iterCount = iterCount +1;   


end
figure;
plot(Fitness);
grid on;

