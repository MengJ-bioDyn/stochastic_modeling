function [ deltaT, deltaX ] = stepTauLeapSingle(x,h,S,tau,epsilon)
% [ deltaT, deltaX ] = stepTauLeap(x,h,S,tau,epsilon)
% deltaT is nan on error

rates = h(x);
combined = sum(rates);


if combined<eps || any(rates<0)
    warning('Rates close to or below zero, stoping simulation');
    deltaT = nan;
    deltaX = nan(size(x));
    return;
end

% Naive approach to choose correct tau
change = 0;
while change <= epsilon*combined
    tau = tau*1.2;
    xprime = x + S*rates*tau;
    change = h(xprime)-rates;
    change = sqrt(sum(change.*change));
end
change=inf;
while change > epsilon*combined
    tau = tau/1.2;
    xprime = x + S*rates*tau;
    change = h(xprime)-rates;
    change = sqrt(sum(change.*change));
end

% Which reactions occur how often
r = poisson(rates*tau);

% Return deltas
deltaX = S*r;
deltaT = tau;


end

