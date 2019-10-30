% Naive arc-length calculation of a curve in R^2 (or C)-plane
%
% mikael.mieskolainen@cern.ch, 2019

function [summ, segs] = arcintegral(x,y)

% Naive arc-length calculation
summ = 0.0;
segs  = zeros(length(x)-1,1);

for i = 1:length(x)-1
    segs(i) = sqrt( (x(i+1)-x(i))^2 + (y(i+1)-y(i))^2 );
    summ = summ + segs(i);
end

end