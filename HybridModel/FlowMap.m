function [Xout,Sout] = FlowMap(X,t,j,S,B,J,noPVInv,measuremedOutputs)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%     dV = B * X(15:15+noPVInv-1);
    dV = zeros(14,1);
    dVS = zeros(28,1);
    if(mod(t,900) == 2)
        dV = B * X(15:15+noPVInv-1);
        if (j>1)
            dS = measuremedOutputs(:,j) - measuremedOutputs(:,j-1);
            dVS = J\dS;
        end
    end
    
    Xout = [X(1:14) + dV + dVS(15:end); X(15:15+noPVInv-1); X(end)+1];
    Sout = S + dVS;
end