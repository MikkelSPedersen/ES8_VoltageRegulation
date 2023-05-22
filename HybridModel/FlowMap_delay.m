function [Xout,Sout] = FlowMap_delay(X,t,j,S,B,J,noPVInv,measuremedOutputs,delay)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%     dV = B * X(15:15+noPVInv-1);
    dV = zeros(14,1);
    dVS = zeros(28,1);
    if(mod(t,900) == 2)
        if (j>1)
            dS = measuremedOutputs(:,j) - measuremedOutputs(:,j-1);
            dVS = J\dS;
        end
    end
    if(mod(t,900) == delay)
        dV = B * X(15:15+noPVInv-1);
    end
    
    Xout = [X(1:14) + dV + dVS(15:end); X(15:15+noPVInv-1); X(end)+1];
    Sout = S + dVS;
end