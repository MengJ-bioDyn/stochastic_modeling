function [ deltaT, deltaX ] = stepLangevinSingle(x,h,S,deltaT)
% [ deltaT, deltaX ] = stepLangevin(x,h,S,deltaT)
% deltaT is nan on error

rates = h(x);
combined = sum(rates);


if combined<eps || any(rates<0)
    warning('Rates close to or below zero, stoping simulation');
    deltaT = nan;
    deltaX = nan(size(x));
    return;
end


% Compute deltaX as the deterministic part ...
deltaX = S*rates*deltaT+...
            ... plus the stochastic part with correct variance:
            S*diag(sqrt(rates))*randn(numel(rates),1)*sqrt(deltaT);

end

