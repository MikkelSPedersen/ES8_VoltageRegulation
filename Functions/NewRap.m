function [Vout,Sout] = NewRap(Sinput)
%NEWRAP Summary of this function goes here
%   Detailed explanation goes here

load(".\DataAndVariables\BaseDefinitions.mat")

% to convert the line admitances from rectangular to polar form
[theta, rho] = cart2pol(real(Y), imag(Y));
% theta = angle
% rho = magnitude

pqBus = 2:15;
nrPQ = length(pqBus);

delta = zeros(nrBus, 1); % voltage phase angle (initialise to slack bus, if possible)
V = ones(nrBus, 1); % voltage magnitude (initialise to slack bus, if possible)

Pcal = zeros(nrBus,1);
Qcal = zeros(nrBus,1);
Psch = Sinput(1:nrPQ);
Qsch = Sinput(nrPQ+1:nrPQ*2);


J1 = zeros(nrBus-1);
J2 = zeros(nrBus-1, nrPQ);
J3 = zeros(nrPQ, nrBus-1);
J4 = zeros(nrPQ);

tolerance=1;

while (tolerance > 1e-3)
% calculation of powers for pq buses
    for i = 1:nrPQ
        m = pqBus(i); % PQ bus indexes, so equation includes only PQ buses
        Pcal(m) = V(m)*sum(rho(m,:).*V.'.*cos(theta(m,:)-delta(m)+delta.'));
        Qcal(m) = -V(m)*sum(rho(m,:).*V.'.*sin(theta(m,:)-delta(m)+delta.'));
    end
    % calculation of mismatches
    dP = Psch - Pcal(2:nrBus);
    dQ = Qsch - Qcal(2:nrBus);
    % mismatch vector
    mismatch = [dP ; dQ];
    
    % calculation of jacobian matrices
    for i = 1:nrBus-1
        m = i+1; % +1 as to exclude slack bus
        for j = 1:nrBus-1
            n = j+1;
            if i == j % diagonal elements for J1
                % subtraction after sum, due to j = i excluded
                YVsum = sum(rho(m,:).*V.'.*sin(theta(m,:)-delta(m)+delta.')) ...
                    -rho(m,m)*V(m)*sin(theta(m,m));
                %deltas ommited in last term, due to delta(m)-delta(m)=0
                J1(i,i) = V(m)*YVsum;
            else % off-diagonal elements for J1
                J1(i,j) = -rho(m,n)*V(m)*V(n)*sin(theta(m,n) ...
                    -delta(m)+delta(n));
            end 
        end
    end
    
    for i = 1:nrBus-1
        m = i+1;
        for j = 1:nrPQ
            n = pqBus(j); % PQ bus indexes, so equation includes only PQ buses
            if i == j % diagonal elements for J2
                YVsum = sum(rho(m,:).*V.'.*cos(theta(m,:)-delta(m)+delta.')) ...
                    -rho(m,m)*V(m)*cos(theta(m,m));
                J2(i,i) = 2*rho(m,m)*V(m)*cos(theta(m,m))+YVsum;
            else % off-diagonal elements for J2
                J2(i,j) = rho(m,n)*V(m)*cos(theta(m,n)-delta(m)+delta(n));
            end 
        end
    end
    
    for i = 1:nrPQ
        m = pqBus(i);
        for j = 1:nrBus-1
            n = j+1;
            if i == j % diagonal elements for J3
                YVsum = sum(rho(m,:).*V.'.*cos(theta(m,:)-delta(m)+delta.')) ...
                    -rho(m,m)*V(m)*cos(theta(m,m));
                J3(i,i) = V(m)*YVsum;
            else % off-diagonal elements for J3
                J3(i,j) = -rho(m,n)*V(m)*V(n)*cos(theta(m,n) ...
                    -delta(m)+delta(n));
            end 
        end
    end
    
    for i = 1:nrPQ
        m = pqBus(i);
        for j = 1:nrPQ
            n = pqBus(j);
            if i == j % diagonal elements for J4
                YVsum = sum(rho(m,:).*V.'.*sin(theta(m,:)-delta(m)+delta.')) ...
                    -rho(m,m)*V(m)*sin(theta(m,m));
                J4(i,i) = -2*rho(m,m)*V(m)*sin(theta(m,m))-YVsum;
            else % off-diagonal elements for J4
                J4(i,j) = -rho(m,n)*V(m)*sin(theta(m,n)-delta(m)+delta(n));
            end 
        end
    end
    
    % initialize jacobian
    J = [J1 J2; J3 J4];

    % calculate correction vector
    correction = J\mismatch;
    
    delta(2:nrBus) = delta(2:nrBus) + correction(1:nrBus-1);
    V(2:nrBus) = V(2:nrBus) + correction(nrBus:end);
    
    tolerance = max(abs(mismatch));
end

Vout = [delta(2:nrBus) ; V(2:nrBus)];
Sout = [Pcal(2:nrBus) ; Qcal(2:nrBus)];
end

