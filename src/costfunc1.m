% Spectrum fit cost function
%
% mikael.mieskolainen@cern.ch, 2019

function [chi2, YHAT, ndf] = costfunc1(param, XVAL)

global YVAL;          % Data y-axis values
if (nargin == 1)
    global XVAL;      % Data x-axis values
elseif (nargin == 2)
    % We got input
end
global FITRANGE;  % Data fit range
global CONSTRAINTS;
global K;         % Number of components
global BWmode;

global PDG_m0;
global PDG_GAMMA0;
global PDG_m0_err;
global PDG_GAMMA0_err;

% Evaluate the fit function
[YHAT, ~, ~, ~] = totfunc(XVAL, param, K, BWmode);

% Enforce fit range
ind = (XVAL >= FITRANGE(1) & XVAL <= FITRANGE(2));

% Chi2
chi2ss = (YVAL(ind) - YHAT(ind)).^2 ./ YVAL(ind);
chi2ss(isinf(chi2ss) | isnan(chi2ss)) = 0;
chi2 = sum(chi2ss);

% Number of Degrees of Freedom, n - p
ndf = sum(chi2ss ~= 0) - length(param);

% -------------------------------------------------------------------------
if (CONSTRAINTS)
% Enforce parameter physicality
[A,m0,GAMMA0,phi, A_BG, gamma_BG, D_BG] = vec2param(param, K);

% Continuum
if     A_BG < 0, chi2 = chi2 + 1e6; end
if gamma_BG < 0, chi2 = chi2 + 1e6; end
if     D_BG < 0, chi2 = chi2 + 1e6; end

% Resonances
for i = 1:length(A)
   if (A(i) < 0)
       chi2 = chi2 + 1e6;
   end
end

% Minimum width (bounded by detector resolution, already)
for i = 1:length(GAMMA0)
   if (GAMMA0(i) < 0.02)
       chi2 = chi2 + 1e6;
   end
end

% How much away from the PDG expectations
NSIGMA = 5;
for i = 1:length(m0)
   if (abs((PDG_m0(i) - m0(i))/PDG_m0(i)) > NSIGMA * PDG_m0_err(i))
       chi2 = chi2 + 1e6;
   end
end
NSIGMA = 10;
for i = 1:length(GAMMA0)
   if (abs((PDG_GAMMA0(i) - GAMMA0(i))/PDG_GAMMA0(i)) > NSIGMA * PDG_GAMMA0_err(i))
       chi2 = chi2 + 1e6;
   end
end
end
% -------------------------------------------------------------------------

end