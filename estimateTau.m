function [ tau ] = estimateTau( x,h,S,epsilon,tauGuess )

tau = tauGuess;

rates = h(x);
combined = sum(rates);

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


end

