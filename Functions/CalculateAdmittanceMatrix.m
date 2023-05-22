function Y = CalculateAdmittanceMatrix(BMva, BVol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

BZ = (BVol*BVol)/(BMva);        % Base impedance 

% Definition of parameters for the transformer
Sb          = BMva;%BMVA is in KVA and Sb must be in W
Vb_mv       = 10e3;%V
Zb          = Vb_mv^2/Sb;%Ohm

Sr_t        = 630e3;         %[kVA]
Vt_mv       = 10e3;          %[V]
Vk_t        = 4;               %[%]
Zt_phase    = Vk_t / 100 * Vt_mv^2 / Sr_t;
trafo_Rpu   = 1 / sqrt(26) * Zt_phase/Zb;
trafo_Xpu   = sqrt(25/26) * Zt_phase/Zb;
trafo       = [trafo_Rpu,trafo_Xpu];

% Define cable parameters, see table 3.1(STC)
Cable1.R         = 0.3339;
Cable1.X         = 0.1163;
Cable1.B         = 328.6337/2*1e-6;
%
Cable2.R         = 0.5163;
Cable2.X         = 0.1213;
Cable2.B         = 278.0747/2*1e-6;
%
Cable3.R         = 1.0327;
Cable3.X         = 0.1264;
Cable3.B         = 252.7951/2*1e-6;

bus2ToBus3        = 0.2522767; % 30111 to 11501
bus3ToBus4        = 0.242335; % 11501 to 11502
bus4ToBus5        = 0.5797393; % 11502 to 11503
bus5ToBus6        = 0.81834586; % 11503 to 11504
bus6ToBus7        = 0.4517369; % 11504 to 11505
bus7ToBus8        = 0.63504136; % 11505 to 11506
bus8ToBus9        = 0.3740655; % 11506 to 11507
bus9ToBus10       = 0.354182; % 11507 to 11510
bus10ToBus11      = 0.1882755; % 11510 to 11511
bus9ToBus12       = 0.80405432; % 11507 to 11512
bus12ToBus13      = 0.70339219; % 11512 to 11513
bus8ToBus14       = 0.5051748; % 11506 to 11508
bus14ToBus15      = 0.1994602; % 11508 to 11509

%         |  From |  To   |  R                         |   X     |                   B/2                    | X'mer   |
%         |  Bus  |  Bus  | OHM                        |  OHM    |                   OHM                    | TAP (a) |
linefeeder = [1      2      trafo(1)                    trafo(2)                    0                         1;
              2      3      bus2ToBus3*Cable1.R         bus2ToBus3*Cable1.X         bus2ToBus3*Cable1.B       1;
              3      4      bus3ToBus4*Cable1.R         bus3ToBus4*Cable1.X         bus3ToBus4*Cable1.B       1;
              4      5      bus4ToBus5*Cable1.R         bus4ToBus5*Cable1.X         bus4ToBus5*Cable1.B       1;
              5      6      bus5ToBus6*Cable2.R         bus5ToBus6*Cable2.X         bus5ToBus6*Cable2.B       1;
              6      7      bus6ToBus7*Cable2.R         bus6ToBus7*Cable2.X         bus6ToBus7*Cable2.B       1;
              7      8      bus7ToBus8*Cable2.R         bus7ToBus8*Cable2.X         bus7ToBus8*Cable2.B       1;
              8      9      bus8ToBus9*Cable2.R         bus8ToBus9*Cable2.X         bus8ToBus9*Cable2.B       1;
              9      10     bus9ToBus10*Cable3.R        bus9ToBus10*Cable3.X        bus9ToBus10*Cable3.B      1;
              10     11     bus10ToBus11*Cable3.R       bus10ToBus11*Cable3.X       bus10ToBus11*Cable3.B     1;
              9      12     bus9ToBus12*Cable3.R        bus9ToBus12*Cable3.X        bus9ToBus12*Cable3.B      1;
              12     13     bus12ToBus13*Cable3.R       bus12ToBus13*Cable3.X       bus12ToBus13*Cable3.B     1;
              8      14     bus8ToBus14*Cable2.R        bus8ToBus14*Cable2.X        bus8ToBus14*Cable2.B      1;
              14     15     bus14ToBus15*Cable2.R       bus14ToBus15*Cable2.X       bus14ToBus15*Cable2.B     1  ];

% Transform to pu
linefeeder(2:end,3:4) = linefeeder(2:end,3:4)./BZ;
linefeeder(2:end,5) = linefeeder(2:end,5)*BZ;

% extracting data from linefeeder
fromBus = linefeeder(:,1);
toBus = linefeeder(:, 2);
R = linefeeder(:,3);
X = linefeeder(:,4);
B = linefeeder(:,5);
a = linefeeder(:,6);

z = R + 1i*X;
y = 1./z;
b = 1i*B;

nBus = max(max(fromBus),max(toBus));
nBranch = length(fromBus);

% Calculate the matrix
Y = complex(zeros(nBus,nBus),zeros(nBus,nBus));
% Calculate the off-diagonal elements
for i=1:nBranch
   Y(fromBus(i), toBus(i)) = - y(i)/a(i);
   Y(toBus(i), fromBus(i)) = Y(fromBus(i), toBus(i)); % Mirror to the other side
end
% Calculate the diagonal elements
for i=1:nBus
    for j=1:nBranch
        if fromBus(j) == i % If the line originates from bus i
            Y(i,i) = Y(i,i) + y(j)/(a(j)^2) + b(j);
        elseif toBus(j) == i  % If the line ends at bus i
            Y(i,i) = Y(i,i) + y(j) + b(j);
        end
    end
end
Y;
end