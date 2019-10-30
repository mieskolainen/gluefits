% Likelihood maximization function
%
% mikael.mieskolainen@cern.ch, 2019

function L = maxlikelihood(x0)

global MDATA;
global K;

% Normalization of the integral
[l,Q] = likelihood(MDATA, x0, K);

logl = log(l/Q);
logl = logl(~isnan(logl) & ~isinf(logl));

% Extended Log-likelihood, where Q is the integral
L = -sum(logl);

end