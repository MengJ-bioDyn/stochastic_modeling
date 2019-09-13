function [ SI, SC, delays, h, endSim ] = DelayedDegradation( Omega )
% X1 is repressor protein, X2 - ClpXP
% X1+X2       --(c)-->  X2

% Rate constants
c = [0.01]';

% Stoichiometry matrices
SI = [ -1
       -1 ];       % update at initiation
SC = [ 0
       1 ];        % update at completion

delays = [0];  % delay between initiation and completion
                   % for every reaction. 0 = no delay
 
% Rates
h = @(x)(repmat(c,[1 size(x,2)])...
          .*[x(1,:).*x(2,:)/Omega]);

% When to end a simulation
% Here: when protein has gone extinct
endSim = @(x)(x(1,:)==0);

end

