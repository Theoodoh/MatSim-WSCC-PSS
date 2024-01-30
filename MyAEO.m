function [BestX,BestF,HisBestFit]=MyAEO(F_index,MaxIt,nPop)


% FunIndex: Index of function.
% MaxIt: The maximum number of iterations.
% PopSize: The size of population.
% PopPos: The position of population.
% PopFit: The fitness of population.
% Dim: The dimensionality of problem.
% C: The consumption factor.
% D: The decomposition factor.
% BestX: The best solution found so far. 
% BestF: The best fitness corresponding to BestX. 
% HisBestFit: History best fitness over iterations. 
% Low: The low bound of search space.
% Up: The up bound of search space.

[Low,Up,Dim]=PssRange(F_index);
PopPos=zeros(nPop,Dim);
PopFit=zeros(1,nPop);
for i=1:nPop
    PopPos(i,:)=rand(1,Dim).*(Up-Low)+Low;
    PopFit(i)=pss_objf(PopPos(i,:),F_index,Dim);
end

[~, indF]=sort(PopFit,'descend');
PopPos=PopPos(indF,:);
PopFit=PopFit(indF);
BestF=PopFit(end);
BestX=PopPos(end,:);

HisBestFit=zeros(MaxIt,1);
Matr=[1,Dim];

for It=1:MaxIt
    r1=rand;
    a=(1-It/MaxIt)*r1;
    xrand=rand(1,Dim).*(Up-Low)+Low;
    newPopPos(1,:)=(1-a)*PopPos(nPop,:)+a*xrand; %equation (1)
    
    u=randn(1,Dim);
    v=randn(1,Dim);
    C=1/2*u./abs(v); %equation (4)
    newPopPos(2,:)=PopPos(2,:)+C.*(PopPos(2,:)-newPopPos(1,:)); %equation (6)
    
    for i=3:nPop        
        u=randn(1,Dim);
        v=randn(1,Dim);
        C=1/2*u./abs(v);
        r=rand;
        if r<1/3
            newPopPos(i,:)=PopPos(i,:)+C.*(PopPos(i,:)-newPopPos(1,:)); %equation (6)
        else
            if  1/3<r && r<2/3
                newPopPos(i,:)=PopPos(i,:)+C.*(PopPos(i,:)- PopPos(randi([2 i-1]),:)); %equation (7)                
            else
                r2=rand;
                newPopPos(i,:)=PopPos(i,:)+C.*(r2.*(PopPos(i,:)- newPopPos(1,:))+(1-r2).*(PopPos(i,:)-PopPos(randi([2 i-1]),:))); %equation (8)                
            end
        end        
    end
    
    for i=1:nPop
        newPopPos(i,:)=SpaceBound(newPopPos(i,:),Up,Low);
        newPopFit(i)=pss_objf(newPopPos(i,:),F_index,Dim);
        if newPopFit(i)<PopFit(i)
            PopFit(i)=newPopFit(i);
            PopPos(i,:)=newPopPos(i,:);
        end
    end
    
    [~, indOne]=min(PopFit);
    for i=1:nPop
        r3=rand;   Ind=round(rand)+1;
        newPopPos(i,:)=PopPos(indOne,:)+3*randn(1,Matr(Ind)).*((r3*randi([1 2])-1)*PopPos(indOne,:)-(2*r3-1)*PopPos(i,:)); %equation (9)
    end
    
    for i=1:nPop
        newPopPos(i,:)=SpaceBound(newPopPos(i,:),Up,Low);
        newPopFit(i)=pss_objf(newPopPos(i,:),F_index,Dim);
        if newPopFit(i)<PopFit(i)
            PopPos(i,:)=newPopPos(i,:);
            PopFit(i)=newPopFit(i);
        end
    end
    
    [~,indF]=sort(PopFit,'descend');    
    PopPos=PopPos(indF,:);
    PopFit=PopFit(indF);
        
    if PopFit(end)<BestF
        BestF=PopFit(end);
        BestX=PopPos(end,:);
    end
    
    HisBestFit(It)=BestF;
    display(['Iteration ' num2str(It) ': HisBestF = ' num2str(HisBestFit(It))]);
end
    