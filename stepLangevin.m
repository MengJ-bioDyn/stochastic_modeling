function [ deltaT, deltaX ] = stepLangevin(x,h,S,deltaT)
% [ deltaT, deltaX ] = stepLangevin(x,h,S,deltaT)
% deltaT is nan on error

rates = h(x);
combined = sum(rates);


errors = combined<eps | any(rates<0,1);
rates(:,errors) = nan;
deltaT(errors) = nan;


% Construct a big version of the stoichiometry matrix
% (repeat blocks of S along the diagonal),
% so it can be applied to multiple states.
global Srepeat;
if any(size(Srepeat)~=[size(x,2)*size(S,1),size(x,2)*size(S,2)])
   if size(Srepeat,1)>size(x,2)*size(S,1) && size(Srepeat,2)>size(x,2)*size(S,2)
       %disp('Shrinking...');
       Srepeat = Srepeat(1:size(x,2)*size(S,1),1:size(x,2)*size(S,2));
   else
       disp('Constructing multi-stoichiometry...');
       Srepeat = spalloc(size(x,2)*size(S,1),size(x,2)*size(S,2),size(x,2)*numel(S));
       for idx=1:size(x,2)
           Srepeat( (1+(idx-1)*size(S,1)):(idx*size(S,1)), ...
                    (1+(idx-1)*size(S,2)):(idx*size(S,2))      ) = S;
       end
   end
end

% Compute deltaX as the deterministic part ...
deltaX = (S*rates).*repmat(deltaT,[size(x,1) 1])+...
         ... plus the stochastic part with correct variance:
         reshape(Srepeat*sparse(1:numel(rates),1:numel(rates),sqrt(rates(:)))*...
                 (randn(numel(rates),1).* ...
                    reshape(repmat(sqrt(deltaT),[size(rates,1) 1]),[numel(rates) 1])),...
                 [size(x,1) size(x,2)]);

end

