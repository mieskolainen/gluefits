% Parameter vector to parameters
%
% mikael.mieskolainen@cern.ch, 2019

function [A,m0,GAMMA0,phi, A_BG,gamma_BG, D_BG] = vec2param(param, K)

% Collect the parameters for Breit-Wigners
A      = zeros(K,1);
m0     = zeros(K,1);
GAMMA0 = zeros(K,1);
phi    = zeros(K,1);

for k = 1:K
   A(k)      = param(1 + (k-1)*4);
   m0(k)     = param(2 + (k-1)*4);
   GAMMA0(k) = param(3 + (k-1)*4);
   phi(k)    = param(4 + (k-1)*4);
end

% Collect the background function parameters
A_BG     = param(end-2);
gamma_BG = param(end-1);
D_BG     = param(end);

end