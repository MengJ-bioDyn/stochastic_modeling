function [ X ] = poisson( lambda )
%[ X ] = poisson( lambda )

if numel(lambda)==0
    X = zeros(size(lambda));
    return;
end

X = zeros(size(lambda));
times = X;
X(isnan(lambda)) = nan;
large = lambda > 20;
if any(large)
    ll = lambda(large);
    X(large) = floor(ll+randn(size(ll)).*sqrt(ll));
    times(large) = 1; % Don't update these anymore
end
while true
    cand = times < 1;
    if all(~cand)
        break;
    end
    l = lambda(cand);
    interval = -log(rand(size(l)))./l;
    times(cand) = times(cand) + interval;
    X(times<1) = X(times<1)+1;
end

end

