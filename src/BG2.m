% Smooth background function (Rice type)
%
% mikael.mieskolainen@cern.ch, 2019

function y = BG2(x, D, gamma, nu)

y = ( (x - D)/gamma^2 ) .* exp( -((x - D).^2 + nu^2) / (2*gamma^2)) .* besseli(0,(x*abs(nu))/gamma^2);


end