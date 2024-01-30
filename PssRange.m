function [Low,Up,Dim]=PssRange(FunIndex)
 
Dim = 10;

switch FunIndex
    case 1
        Low = [0.001,0.001,0.001,0.02,0.001,0.02,0.001,0.02,0.001,0.02];
        Up = [50,50,1,1,1,1,1,1,1,1];
    otherwise
end