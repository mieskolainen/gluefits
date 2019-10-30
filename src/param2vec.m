% Parameters to a parameter vector
% 
% mikael.mieskolainen@cern.ch, 2019

function param = param2vec(A,m0,GAMMA0,phi, A_BG,gamma_BG,D_BG)

K = length(A);
param = zeros(4*K + 2,1);

% Collect the parameters for Breit-Wigners
for k = 1:K
    param(1 + (k-1)*4) = A(k);
    param(2 + (k-1)*4) = m0(k);
    param(3 + (k-1)*4) = GAMMA0(k);
    param(4 + (k-1)*4) = phi(k);
end

% Collect the background function parameters
param(end-2) = A_BG;
param(end-1) = gamma_BG;
param(end)   = D_BG;

end