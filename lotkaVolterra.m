function [ S, h, endSim ] = lotkaVolterra( Omega )
% X1 is prey
% X2 is predator
% X1       --(c1)-->  2X1
% X1 + X2  --(c2)-->  2X2
% X2       --(c3)-->    0

% Rate constants
c = [0.5 0.0025 0.3]';

% Stoichiometry matrix
S = [1  -1  0
     0   1  -1];

 
% Rates
h = @(x)(repmat(c,[1 size(x,2)])...
          .*[x(1,:)
             x(1,:).*x(2,:)/Omega
             x(2,:)]);

% When to end a simulation
% Here: when predator has gone extinct
endSim = @(x)(x(2,:)==0);

end

