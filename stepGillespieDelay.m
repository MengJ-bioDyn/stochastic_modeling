function [ deltaT, deltaX ] = stepGillespieDelay(x,t,h,SI,SC,delays)
% [ deltaT, deltaX ] = stepGillespie(x,h,S)
% deltaT is nan on error

global Q R nq;

rates = h(x);
combined = sum(rates);
csrates = cumsum(rates)/combined;

if (combined<eps && nq==0) || any(rates<0)
    warning('Rates below zero, stoping simulation');
    deltaT = nan;
    deltaX = nan(size(x));
    return; 
end

% Find time step
deltaT = -log(rand)/(combined);

% Compare t+deltaT with the queue
if nq == 0 || t+deltaT < min(Q(1:nq))
    % Find reaction
    reaction = find(rand<=csrates,1,'first'); 

    % Check reaction type
    if delays(reaction) == 0  %ND - No delay
    % Choose right deltaX from stoichiometry matrix
        deltaX = SC(:,reaction)+SI(:,reaction);
    else % There is a delay between initiation and completion, CD or ICD
        deltaX=SI(:,reaction);
        nq=nq+1;
        Q(nq)=t+deltaT+delays(reaction);
        R(nq)=reaction;
        [Q,IX]=sort(Q,'descend');
        R = R(IX);
    end
    
else
    deltaT=Q(nq)-t;
    Q(nq) = -1;
    deltaX=SC(:,R(nq));
    R(nq) = 0;
    nq=nq-1;
end  

end

