function [Xout, Sout, k, Bout, Jout, Uout] = JumpMap(X, S, Kin, Uin, Hu, Hp, Q, R, rho, uMax, measuremedOutputs, ConsumptionEst, noPVInv, PVIndex, ConsumptionEstData)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    
    % Linearise
    [Vout, SLin, Sensitivity, J] = NewRapV2(S);
    B = Sensitivity(:,PVIndex);
    Bout = B;
    Jout = J;
    A = eye(14);
    ref = ones(14,1) - Vout(15:end);
    Vlin = Vout(15:end);

    % Calculate next predicted consumption change
    dVest = zeros(28,1);
    k = Kin+1;
    if (k <= 96)
        dSest = (ConsumptionEstData(:,k) - ConsumptionEstData(:,k-1));
        dVest = J\dSest;
    end

    %MPC
    uopt = sdpvar(noPVInv, Hu);
    control = sdpvar(noPVInv, Hp);
    duopt = sdpvar(noPVInv, Hp);
    slack = sdpvar(2,Hp);
    constraints = [];
    objective = 0;

    xopt = X(1:14) - Vlin;
    % This loop iterates over the entire prediction horizon
    for j = 1:Hp
        if (j <= Hu) % Checks if we are below our control horizon
            control(:,j) = uopt(:,j);
        else
            control(:,j) = uopt(:,end);
        end
        if (j == 1)
            duopt(:,j) = control(:,j) - Uin;
        else
            duopt(:,j) = control(:,j) - control(:,j-1);
        end
        if (j == 1 && ConsumptionEst == true)
            xopt = A * xopt + B * duopt(:,j) + dVest(15:end);
        else
            xopt = A * xopt + B * duopt(:,j);
        end
        objective = objective + (xopt-ref).' * Q *(xopt-ref) + duopt(:,j).'* R *duopt(:,j) +  rho*(slack(1,j) + slack(2,j));
        constraints = [constraints,  0.95 + slack(1,j) <= Vlin+xopt <= 1.05+ slack(2,j)];
    end
    ops = sdpsettings('solver','ipopt');
    ops.verbose = 0; % makes yalmip quiet
    % solving the problem
    optimize([constraints, slack >= 0, abs(uopt) <= uMax],objective, ops);
    
%     Xout(15:17) = uopt(:,1);
    %duPlot(:,i) = duopt(:,1);
               
    %update sample
    du = double(duopt(:,1));

    Sout = S;
    Uout = uopt(:,1);
    Xout = [X(1:14); du; 1];
end