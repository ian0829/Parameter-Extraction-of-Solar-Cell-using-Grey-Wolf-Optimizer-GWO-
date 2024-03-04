clear
clc
tic
% typical GA
% evaluated with test functions (maxF3)

%% ------------------------------------------
%*************************************%
%---Scenario setup ---%
%*************************************%
% HW:ø�X��{���ϧ�
% Hint: ���O meshgrid()/surf()

xx = -1:0.1:3
yy = -3:0.1:3
[x,y] = meshgrid(xx, yy);

term1 =  20.*exp((-0.2).*(x.^2+y.^2));
term2 =  exp(cos(2.*pi.*x)+cos(2.*pi.*y));
z = term1 + term2 ;

figure;
surf(x, y, z);

%---------------------------------------
% ----------- GA Parameters ------------
%----------------------------------------
Generation = 100;   % 100 �@�N
NC = 100;           % 100��chrmosome
XbitN = 12;         % 12 bits�N��x  ( 1 < x < 3) 
YbitN = 13;         % 13 bits�N��y  (-3 < y < 3) 
ChromoBit = XbitN + YbitN ; % �C��chromosome��bits��
popu = zeros( NC,ChromoBit ); % population ( �x�}�j�p 100x25 )
Fitness = zeros( NC, 1 ); % �إߤ@�x�}��m�A�s��
best = zeros(Generation,1);  % �إߤ@�x�}��m�ӥ@�N���u
bestR = zeros(Generation,2); % �إߤ@�x�}��m�ӥ@�N���u��x,y��
most_best = zeros(Generation,1); % �إߤ@�x�}��m�ðê�
ratio = zeros( NC, 1 ); % �إߤ@�x�}��m�Q����v

%% --------------------------------------------------
%*************************************%
%--- �� �l ��  Initialization ---%
%*************************************%
% HW: �إߤj�p��100x25����]�x�}(100 chrmosome, �C��chromosome��25��Bits)
% Hint: 2�h for �j��/ ���Ͷü�[0,1]�è����
%

for i = 1:100
    for j = 1:25
        popu(i, j) = round(rand);
    end
end

%% --------------------------------------------------
for genCounter = 1 : Generation
 %% --------------------------------------------------   
    %*************************************%
    % ---Evalution--- %
    %*************************************%
 % HW: �إߤj�p��100x25����]�x�}(100 chrmosome, �C��chromosome��25��Bits)
 % Hint: step 1 : �Nx�q�G�i���ন�Q�i��
 % Hint: step 2 : �u���ഫ�A�����ܯu��y�нd��
 % Hint: step 3 : �o��u��x,y�ƭ� �ña�J������{��
    
 for i = 1:NC
      p = bi2de( popu(i, 1:12) );
      P(i) = -1 + (p./(2^12-1)).*4
      q = bi2de(popu(i, 13:25));
      Q(i) = -3 + (q./(2^13-1)).*6
      R(i, 1) = P(i);
      R(i, 2) = Q(i);
      
      Fitness(i, 1) = maxF3(R(i,:));
 end 
  
 %% --------------------------------------------------   
    % �������@�N�̦n����(best) �H�� �����L�h�@�N��{�b�̦n����(most_best)
     best( genCounter, 1 )= max( Fitness );
    if genCounter == 1
       most_best(1,1)=best(1,1);
    else      
        if best( genCounter,1 ) > most_best( genCounter-1,1 )
            most_best( genCounter,1) = best( genCounter ,1 );
        else 
            most_best( genCounter,1 ) = most_best( genCounter-1 ,1 );
        end
    end 
 %% --------------------------------------------------   
    %****************************************%   
    % ---Selcet Parents[ Roulette wheel ]--- %
    %****************************************%
     % HW: Roulette wheel ��ܿ˥N
    % �]���D�̤j�ȡAFitness�ȶV�j�A�ҥe��V�j�N��V���u��
    % �]����"����"�������CFitness �M �⤤���v������
      
    % step 1 : ���N Fitness �`�M�[�_�ӡA�@������  
    SumFitness = sum(Fitness);
    
    % step 2: ���l���ӧO��Fitness�A�Φ��⤤���ratio  
    ratio(i, 1) = (Fitness(i, 1)) / SumFitness
    
    % step 3 : �}�l�e���
    label(1, 1) = ratio(1, 1)
        for i = 2:100
            label(i, 1) = label (i-1, 1) + ratio (i, 1);
        end
 
   % �}�l����L�i����
   % ��� ���ʧ@�Ѷü�rand(1)�Ӽ���
   for k = 1 : 2 : 100  
        
      [m1,index1] = min(abs(rand(1)-label)); % �P���@�Ө�׳̱���
      chromoFather = popu(index1,1:25);  % �qRoulette wheel�H����X�˥N(����)
      [m2,index2] = min(abs(rand(1)-label));  % �P���@�Ө�׳̱���
      chromoMother = popu(index2,1:25);  % �qRoulette wheel�H����X�˥N(����)

      
      %******************************************%
      % -------------Recombination-------------- % 
      %******************************************%
      % one-point crossover and one-bit mutation %

      crossoverRate = 0.6;
      cross_counter = 80;
      not_cross_counter = 20;
      

       %----- crossover------ %
        if cross_counter > 0 && not_cross_counter > 0
        if rand(1) < crossoverRate
             
            cross_counter = cross_counter-2;

             cross_location = ceil(rand(1)*25); %�H�����ͭn�o��crossover���_�I
             kid1( 1,1:cross_location )  = chromoFather( 1, 1: cross_location );
             kid1( 1,cross_location+1:25 ) = chromoMother( 1, cross_location+1:25 );
             kid2( 1,1:cross_location )  = chromoMother( 1, 1: cross_location );
             kid2( 1,cross_location+1:25 ) = chromoFather( 1, cross_location+1:25 );

        else
 

        popu( k, 1:25 )   = kid1;  %crossover��զ����l�N
        popu( k+1, 1:25 ) = kid2; %crossover��զ����l�N

        end
        end
   end
 %% --------------------------------------------------  
  %----- mutation-----%

     mutationRate = 0.05;  
     muta_counter = 0.05*100*25; %�n���ܪ�bits��
 
     % �H���q2500(=100*25)��bits�̿�125��bits�Ӭ���
     for i = 1 : muta_counter
         mu_colum=ceil(rand(1)*25);
         mu_row=ceil(rand(1)*100);
         popu(mu_row,mu_colum)=~popu(mu_row,mu_colum);
     end
 
% A=popu;

end

%% 
% HW: �̫�e�X���Ħ��u
% most_best
% best
toc
figure;

plot(most_best, 'r')
hold on;
polt(best, 'b')
legend('Best record so far', 'Best in population', 'Location', 'sotheast');
xlabel('Generation');
ylabel('Fitness')
title('Convergence History')










     
     
    
  
  
  
  
  
  
  
  
        