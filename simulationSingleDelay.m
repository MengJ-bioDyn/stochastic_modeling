clear global;
global Q R nq; 

% % Degrade-and-Fire oscillator
% % Remove comment marker from the following block
% % to simulate this model
% %///////////////////////////////
% % Volume
% Omega = 1;
% 
% % Get reactions
% [ SI, SC, delays, h, endSim ] = DegradeAndFire( Omega );
% 
% % Initial condition
% x0 = round(Omega*5);
% 
% % Simulation time
% Tmax = 50; 
% %///////////////////////////////



% Delayed degradation
% Remove comment marker from the following block
% to simulate this model
%///////////////////////////////
% Volume
Omega = 1;

% Get reactions
[ SI, SC, delays, h, endSim ] = DelayedDegradation( Omega );

% Initial condition
x0 = round(Omega*[500 10]');   % initial protein and ClpXP

% Simulation time
Tmax = 500; 
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
Q(1:1000)=-100;  % array of delayed times (in descending order)
R(1:1000)=0;  % array of delayed reactions
nq=0;         % number of queued reactions


while t<Tmax
    
    % Direct Gillespie
    [deltaT, deltaX] = stepGillespieDelay(x,t,h,SI,SC,delays);
    
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
plot(Treg,X/Omega,'-');
xlabel('t');
ylabel('concentrations');


