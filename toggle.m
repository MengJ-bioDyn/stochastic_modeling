function [ S, h, endSim ] = toggle( Omega )
% 0  --(k1/(1+([X2]/Ki1)^p1)--> X1
% X1 --(k2)-->              0
% 0  --(k3/(1+([X1]/Ki2)^p2)--> X1
% X2 --(k4)-->              0

% Rate prefactors
k1 = 1;
k3 = 1;

k2 = 0.1;
k4 = k2;

Ki1 = 3;
Ki2 = Ki1;

p1=2;
p2=p1;

% Rate constants 
c = @(x)([k1./(1+(x(2,:)/(Ki1*Omega)).^p1)
          k2*ones(1,size(x,2))
          k3./(1+(x(1,:)/(Ki2*Omega)).^p2)
          k4*ones(1,size(x,2))]);

% Stoichiometry matrix
S = [ 1 -1  0  0
      0  0  1 -1 ];

% Rates
h = @(x)(c(x).*[ Omega*ones(1,size(x,2))
                 x(1,:)
                 Omega*ones(1,size(x,2))
                 x(2,:) ]);

% Never end the simulation
endSim = @(x)(false(1,size(x,2)));


end

