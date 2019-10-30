% Smooth background function
%
% mikael.mieskolainen@cern.ch, 2019

function y = BG1(x, gamma, D)

y = ( (x - D) /gamma ) .* exp( -(x - D) / gamma);
y = y.^(1/2);

end