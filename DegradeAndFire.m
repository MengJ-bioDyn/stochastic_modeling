function [ SI, SC, delays, h, endSim ] = DegradeAndFire( Omega )
% X is repressor protein
% 0       --(c1/(1+[Xt]^2)-->  X
% X       --(c2/(1+X))-->    0

% Rate constants
c = [180 20.]';

% Stoichiometry matrix in
SI = [0  -1];     % update at initiation
SC = [1  0];      % update at completion
delays = [1 0];   % delay between initiation and completion
                  % for every reaction. 0 = no delay
 
% Rates
h = @(x)(repmat(c,[1 size(x,2)])...
          .*[1./(1+(x(1,:)/Omega))^2
             x(1,:)./(1+x(1,:))]);

% When to end a simulation
% Here: when predator has gone extinct
endSim = @(x)(x(2,:)==0);

end

