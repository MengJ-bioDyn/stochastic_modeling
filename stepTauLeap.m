function [ deltaT, deltaX ] = stepTauLeap(x,h,S,tau,epsilon,maxTau)
% [ deltaT, deltaX ] = stepTauLeap(x,h,S,tau,epsilon)
% deltaT is nan on error
rates = h(x);
combined = sum(rates);

if ~exist('maxTau','var')
    maxTau = inf;
end
    
errors = combined<eps | any(rates<0,1);
rates(:,errors) = nan;
tau(errors) = nan;

% Naive approach to choose correct tau
needupdate = ~errors;
while any(needupdate)
    tau(needupdate) = tau(needupdate)*1.7;
    xprime = x(:,needupdate) + repmat(tau(needupdate),[size(x,1) 1]).*(S*rates(:,needupdate));
    change = h(xprime)-rates(:,needupdate);
    change = sqrt(sum(change.*change,1));
    
    oldneedupdate = needupdate;
    needupdate = false(1,size(x,2));
    needupdate(oldneedupdate) = change <= epsilon*combined(oldneedupdate) ...
                                & tau(oldneedupdate) < maxTau;
end
change=inf;
needupdate = ~errors;
while any(needupdate)
    tau(needupdate) = tau(needupdate)/1.7;
    
    xprime = x(:,needupdate) + repmat(tau(needupdate),[size(x,1) 1]).*(S*rates(:,needupdate));
    change = h(xprime)-rates(:,needupdate);
    change = sqrt(sum(change.*change));
    
    oldneedupdate = needupdate;
    needupdate = false(1,size(x,2));
    needupdate(oldneedupdate) = (change > epsilon*combined(oldneedupdate)) ...
                                 & (tau(oldneedupdate) > 2./combined(oldneedupdate));
end

deltaX = nan(size(x));


% If time tau step is too small, do a Gillespie step
gillespie = (tau < 2./combined) & (~errors);
if sum(gillespie) > 0
    [tau(gillespie),deltaX(:,gillespie)] = stepGillespie(x(:,gillespie),h,S);
end

if sum(~gillespie) > 0
    % Which reactions occur how often
    r = poisson(rates(:,~gillespie).*(repmat(tau(~gillespie),[size(rates,1) 1])));
    deltaX(:,~gillespie) = S*r;
end

% Return deltas
deltaT = tau;


end

