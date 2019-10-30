% Smooth background function
%
% mikael.mieskolainen@cern.ch, 2019

function f = BG3(x, D, alpha, beta)

f = beta^alpha*(x - D).^(alpha-1).*exp(-(x - D)*beta) / gamma(alpha);

end