

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


%% Simulation (multiple)
% Choose the simulation method by commenting in/out the corresponding
% lines in the loop

reset(RandStream.getGlobalStream);

% Number of trajectories to simulate
numReals = 1000;

t = repmat(0,[1 numReals]);
x = repmat(x0,[1 numReals]);

idx = 0;

tsampleIdx = ones(1,numReals);
deltaSample = 0.1;
X = nan(size(x,1),ceil(Tmax/deltaSample+1),size(x,2));
Treg = (0:size(X,2)-1)*deltaSample;
T = nan(1,size(X,2),size(x,2));
X(:,1,:) = x;
T(1,1,:) = t;

% for tau leap:
tau = 0.01*ones(1,size(x,2));  %initial tau
epsilon = 0.01;                %accuracy constraint

% for Langevin and deterministic:
deltaTFixed = 0.1;
clearvars -global Srepeat;
    

while any(t<Tmax)
    updateThese = t<Tmax;
    
    if mod(idx,100)==0
        disp([num2str(sum(updateThese)) ' to update, median(t) = ' num2str(median(t(updateThese)))]);
    end
    idx = idx + 1;
    
    % Direct Gillespie
    [deltaT, deltaX] = stepGillespie(x(:,updateThese),h,S);
    
    % Tau leap
    %[deltaT, deltaX] = stepTauLeap(x(:,updateThese),h,S,tau(updateThese),epsilon);
    %tau(updateThese) = deltaT;
    
    % Langevin
    %[deltaT, deltaX] = stepLangevin(x(:,updateThese),h,S,deltaTFixed*ones(1,sum(updateThese)));
    
    % Deterministic
    %deltaX = deltaTFixed*S*h(x);
    %deltaT = repmat(deltaTFixed,[1 size(x,2)]);
    

    
    % Update state and time
    lastx = x;
    lastt = t;
    
    x(:,updateThese) = x(:,updateThese) + deltaX;
    t(updateThese) = t(updateThese) + deltaT;
    
    % For some algorithms, molecule numbers can dip below zero
    x(x<0) = 0;
    
    while(true)
        saveThese = t>(1+eps)*tsampleIdx*deltaSample & tsampleIdx+1<=size(X,2);
        if(~any(saveThese))
            break;
        end
        tsampleIdx(saveThese) = tsampleIdx(saveThese)+1;
        for saveThis = find(saveThese)
            X(:,tsampleIdx(saveThis),saveThis) = lastx(:,saveThis);
            T(1,tsampleIdx(saveThis),saveThis) = lastt(saveThis);
        end
    end
    
    % End some simulations
    endThese = endSim(x);
    x(:,endThese) = nan;
    t(endThese)=nan;
    
end

%figure;

% Molecule Numbers
%plot(Treg,squeeze(X(1,:,:)));
%xlabel('t');
%ylabel('molecule numbers');

% Concentrations
%plot(Treg,squeeze(X(1,:,:))/Omega);
%xlabel('t');
%ylabel('concentrations');

%% Plot statistics
figure;
colors = {'b','g','r'};

scaleFactor = 1;

% Comment out for molecule numbers
scaleFactor = 1/Omega;

for idx=1:size(X,1)
    trace = nanmean(X(idx,:,:),3);
    lower = prctile(X(idx,:,:),2.5,3);
    upper = prctile(X(idx,:,:),97.5,3);
    
    
    hp=plot(Treg,trace*scaleFactor,'LineWidth',2,'Color',colors{idx});
    hold off;
    hold on;
    plot(Treg,[lower; upper]*scaleFactor,'--','Color',get(hp,'Color'));
    hold all;
end
hold off;

xlabel('t');
ylabel('mean and percentiles');