% Likelihood function
%
% mikael.mieskolainen@cern.ch, 2019

function [l,Q] = likelihood(x, param, K)

global BWmode;

% Normalization integral
m_min = 0.28;
m_max = 100;
Q = integral(@(m)totfunc(m,param,K,BWmode), m_min, m_max);
l = totfunc(x,param,K,BWmode);

end