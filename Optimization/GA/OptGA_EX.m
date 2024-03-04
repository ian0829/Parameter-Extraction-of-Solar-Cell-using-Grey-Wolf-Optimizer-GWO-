clear
clc
tic
% typical GA
% evaluated with test functions (maxF3)

%% ------------------------------------------
%*************************************%
%---Scenario setup ---%
%*************************************%
% HW:繪出方程式圖形
% Hint: 指令 meshgrid()/surf()

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
Generation = 100;   % 100 世代
NC = 100;           % 100個chrmosome
XbitN = 12;         % 12 bits代表x  ( 1 < x < 3) 
YbitN = 13;         % 13 bits代表y  (-3 < y < 3) 
ChromoBit = XbitN + YbitN ; % 每個chromosome的bits數
popu = zeros( NC,ChromoBit ); % population ( 矩陣大小 100x25 )
Fitness = zeros( NC, 1 ); % 建立一矩陣放置適存值
best = zeros(Generation,1);  % 建立一矩陣放置該世代最優
bestR = zeros(Generation,2); % 建立一矩陣放置該世代最優的x,y值
most_best = zeros(Generation,1); % 建立一矩陣放置衛冕者
ratio = zeros( NC, 1 ); % 建立一矩陣放置被選機率

%% --------------------------------------------------
%*************************************%
%--- 初 始 化  Initialization ---%
%*************************************%
% HW: 建立大小為100x25的基因矩陣(100 chrmosome, 每個chromosome有25個Bits)
% Hint: 2層 for 迴圈/ 產生亂數[0,1]並取整數
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
 % HW: 建立大小為100x25的基因矩陣(100 chrmosome, 每個chromosome有25個Bits)
 % Hint: step 1 : 將x從二進制轉成十進制
 % Hint: step 2 : 線性轉換，對應至真實座標範圍
 % Hint: step 3 : 得到真實x,y數值 並帶入評價方程式
    
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
    % 紀錄本世代最好的值(best) 以及 紀錄過去世代到現在最好的值(most_best)
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
     % HW: Roulette wheel 選擇親代
    % 因為求最大值，Fitness值越大，所占比越大代表越具優勢
    % 因此為"正比"的概念。Fitness 和 抽中機率成正比
      
    % step 1 : 先將 Fitness 總和加起來，作為分母  
    SumFitness = sum(Fitness);
    
    % step 2: 分子為個別的Fitness，形成抽中比例ratio  
    ratio(i, 1) = (Fitness(i, 1)) / SumFitness
    
    % step 3 : 開始畫刻度
    label(1, 1) = ratio(1, 1)
        for i = 2:100
            label(i, 1) = label (i-1, 1) + ratio (i, 1);
        end
 
   % 開始轉輪盤進行選擇
   % 轉動 的動作由亂數rand(1)來模擬
   for k = 1 : 2 : 100  
        
      [m1,index1] = min(abs(rand(1)-label)); % 與哪一個刻度最接近
      chromoFather = popu(index1,1:25);  % 從Roulette wheel隨機選出親代(爸爸)
      [m2,index2] = min(abs(rand(1)-label));  % 與哪一個刻度最接近
      chromoMother = popu(index2,1:25);  % 從Roulette wheel隨機選出親代(媽媽)

      
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

             cross_location = ceil(rand(1)*25); %隨機產生要發生crossover的斷點
             kid1( 1,1:cross_location )  = chromoFather( 1, 1: cross_location );
             kid1( 1,cross_location+1:25 ) = chromoMother( 1, cross_location+1:25 );
             kid2( 1,1:cross_location )  = chromoMother( 1, 1: cross_location );
             kid2( 1,cross_location+1:25 ) = chromoFather( 1, cross_location+1:25 );

        else
 

        popu( k, 1:25 )   = kid1;  %crossover後組成的子代
        popu( k+1, 1:25 ) = kid2; %crossover後組成的子代

        end
        end
   end
 %% --------------------------------------------------  
  %----- mutation-----%

     mutationRate = 0.05;  
     muta_counter = 0.05*100*25; %要突變的bits數
 
     % 隨機從2500(=100*25)個bits裡選125個bits來突變
     for i = 1 : muta_counter
         mu_colum=ceil(rand(1)*25);
         mu_row=ceil(rand(1)*100);
         popu(mu_row,mu_colum)=~popu(mu_row,mu_colum);
     end
 
% A=popu;

end

%% 
% HW: 最後畫出收斂曲線
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










     
     
    
  
  
  
  
  
  
  
  
        