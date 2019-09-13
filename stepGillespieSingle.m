function [ deltaT, deltaX ] = stepGillespieSingle(x,h,S)
% [ deltaT, deltaX ] = stepGillespie(x,h,S)
% deltaT is nan on error

% rates = h(x);
rates = h;
combined = sum(rates);
csrates = cumsum(rates)/combined;

if combined<eps || any(rates<0)
    warning('Rates close to or below zero, stoping simulation');
    deltaT = nan;
    deltaX = nan(size(x));
    return;
end

% Find time step
deltaT = -log(rand)/combined;

% Find reaction
reaction = find(rand<=csrates,1,'first');

% Chose right deltaX from stoichiometry matrix
deltaX = S(:,reaction);

end

