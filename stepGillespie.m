function [ deltaT, deltaX ] = stepGillespie_v2(x,h,S)
% [ deltaT, deltaX ] = stepGillespie(x,h,S)
% deltaT is nan on error

rates = h(x);
combined = sum(rates,1);
csrates = cumsum(rates,1)./repmat(combined,[size(rates,1) 1]);

%warning('Rates close to or below zero, stoping simulation');
errors = combined<eps | any(rates<0,1);
deltaT(errors) = nan(1,sum(errors));
deltaX(:,errors) = nan(size(x,1),sum(errors));

% Find time step
deltaT(~errors) = -log(rand(1,sum(~errors)))./combined(~errors);

% Find reaction
chooseReactMat = repmat(rand(1,sum(~errors)),[size(csrates,1) 1])<=csrates(:,~errors);
indcs = repmat((1:size(chooseReactMat,1))',[1 size(chooseReactMat,2)]);
indcs(~chooseReactMat) = inf;

reaction = min(indcs,[],1);

% Chose right deltaX from stoichiometry matrix
deltaX(:,~errors) = S(:,reaction);

end

