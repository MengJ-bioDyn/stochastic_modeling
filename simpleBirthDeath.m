function [ S, h, endSim ] = simpleBirthDeath( Omega )
% R is mRNA
% P is protein
% 0   --(alphaR)-->  R
% R   --(betaR )-->  0
% R   --(alphaP)-->  R + P
% P   --(betaP )-->  0

alphaR = 1;
alphaP = 1;
betaR = 0.5;
betaP = 0.5;

% Rate constants
c = [alphaR betaR alphaP betaP]';

% Stoichiometry matrix
S = [1  -1  0   0
     0   0  1  -1];

 
% Rates
h = @(x)(repmat(c,[1 size(x,2)])...
          .*[ones(1,size(x,2))*Omega
             x(1,:)
             x(1,:)
             x(2,:)]);

% When to end a simulation
% Here: never
endSim = @(x)(false(1,size(x,2)));

end

