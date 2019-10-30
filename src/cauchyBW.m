% Non-relativistic Breit-Wigner (Cauchy distribution) as a function of m
%
% mikael.mieskolainen@cern.ch, 2019

function f = cauchyBW(m, m0, GAMMA0)

N = pi*0.5*GAMMA0;
f = (1/N)*(0.5*GAMMA0)^2./((m - m0).^2 + (0.5*GAMMA0)^2);

end
