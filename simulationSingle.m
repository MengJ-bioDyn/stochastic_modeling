

% Lotka-Volterra
% Remove comment marker from the following block
% to simulate this model
%///////////////////////////////
% % Volume
% Omega = 1;
% 
% % Get reactions
% [ S, h, endSim ] = lotkaVolterra( Omega );
% 
% % Initial condition
% x0 = round(Omega*[100 100]');
% 
% % Simulation time
% Tmax = 50;
%///////////////////////////////


% Simple birth-death process
% Remove comment marker from the following block
% to simulate this model
%///////////////////////////////
% % Volume
% Omega = 100;
% 
% % Get reactions
% [ S, h, endSim ] = simpleBirthDeath( Omega );
% 
% % Initial condition
% x0 = round(Omega*[0 0]');
% 
% % Simulation time
% Tmax = 50;
%///////////////////////////////


% Goodwin oscillator
% Remove comment marker from the following block
% to simulate this model
%///////////////////////////////
% Volume
Omega = 50;

% Get reactions
[ S, h, endSim ] = goodwin( Omega );

% Initial condition
x0 = round(Omega*[0 0 0]');

% Simulation time
Tmax = 500;
%///////////////////////////////


% Toggle switch
% Remove comment marker from the following block
% to simulate this model
%///////////////////////////////
% % Volume
% Omega = 3;
% 
% % Get reactions
% [ S, h, endSim ] = toggle( Omega );
% 
% % Initial condition
% x0 = round(Omega*[0 5]');
% 
% % Simulation time
% Tmax = 5000;
%///////////////////////////////


%% Simulation
% Choose the simulation method by commenting in/out the corresponding
% lines in the loop.
% In the loop the "...Single" step functions are used, but it would work
% the same without this suffix (just the functions look more complicated,
% because they can handle any number of realizations).

reset(RandStream.getGlobalStream);

t = 0;
x = x0;

idx = 0;

tsampleIdx = 1;
deltaSample = 0.1;
X = nan(numel(x),ceil(Tmax/deltaSample+1));
Treg = (0:size(X,2)-1)*deltaSample;
T = nan(1,size(X,2));
X(:,1) = x;
T(1) = t;

% for tau leap:
tau = 0.01;      %initial tau
epsilon = 0.01;  %accuracy constraint

% for Langevin and deterministic:
deltaTFixed = 0.1;


while t<Tmax
    
    % Direct Gillespie
    [deltaT, deltaX] = stepGillespieSingle(x,h,S);
    
    % Tau leap
    %[deltaT, deltaX] = stepTauLeapSingle(x,h,S,tau,epsilon);
    %tau = deltaT;
    
    % Langevin
    %[deltaT, deltaX] = stepLangevinSingle(x,h,S,deltaTFixed);
    
    % Deterministic
    %deltaX = deltaTFixed*S*h(x);
    %deltaT = repmat(deltaTFixed,[1 size(x,2)]);
    
    % Check if step succeeded
    if isnan(deltaT)
        break;
    end
        
    % Update state and time
    lastx = x;
    lastt = t;
    x = x + deltaX;
    t = t + deltaT;
    
    % For some algorithms, molecule numbers can dip below zero
    x(x<0) = 0;
    
    idx = idx + 1;
    if mod(idx,10000)==0
        t
    end
    
    
    while(t>=(1+eps)*tsampleIdx*deltaSample && tsampleIdx+1<=size(X,2))
        tsampleIdx = tsampleIdx+1;
        X(:,tsampleIdx) = lastx;
        T(tsampleIdx) = lastt;
    end
end

figure;

% Molecule numbers
%plot(Treg,X);
%xlabel('t');
%ylabel('molecule numbers');

% Concentrations
plot(Treg,X/Omega);
xlabel('t');
ylabel('concentrations');

