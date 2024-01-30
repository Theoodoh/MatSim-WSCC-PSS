%--------------------------------------------------------------------------
% Artificial ecosystem-based optimization (AEO)
% AEO code v1.1.
% Developed in MATLAB R2011b
% The code is based on the following papers.
% Zhao, W., Wang, L., & Zhang, Z. (2020). Artificial ecosystem-based 
% optimization: a novel nature-inspired meta-heuristic algorithm. Neural 
% Computing and Applications, 32(13), 9383-9425.
% -------------------------------------------------------------------------

clc;
clear;
 
MaxIteration=100; 
PopSize=100;
  FunIndex=1;
[BestX,BestF,HisBestF]=MyAEO(FunIndex,MaxIteration,PopSize);


display(['F_index=', num2str(FunIndex)]);
display(['The best fitness is: ', num2str(BestF)]);
display(['The best solution is: ', num2str(BestX)]);


figure;
%plot(HisBestF,'r','LineWidth',2);
semilogy(HisBestF,'r','LineWidth',2);
xlabel('Iteration');
ylabel('Fitness');
title(['F',num2str(FunIndex)]);
grid on;







