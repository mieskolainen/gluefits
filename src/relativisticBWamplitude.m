% Relativistic Breit-Wigner as a function of m (not m^2)
%
% mikael.mieskolainen@cern.ch, 2019

function f = relativisticBWamplitude(m, m0, GAMMA0, m_daughter)

% Momentum dependent width
if (m_daughter > 0)
    GAMMA = GAMMA0 * (m0./m) .* ((m.^2 - 4*m_daughter^2)./(m0^2 - 4*m_daughter^2)).^(3/2);
else
    GAMMA = GAMMA0;
end

% sqrt(2m) comes from |f(m^2)|^2 dm^2 -> |f(m)|^2 dm Jacobian
f = sqrt(2*m) ./ (m.^2 - m0^2 + 1i*m0*GAMMA);

end
