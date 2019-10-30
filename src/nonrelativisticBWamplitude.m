% Non-relativistic Breit-Wigner amplitude as a function of m
%
% http://pdg.lbl.gov/2016/reviews/rpp2016-rev-resonances.pdf
%
% mikael.mieskolainen@cern.ch, 2019

function f = nonrelativisticBWamplitude(m, m0, GAMMA0, m_daughter)

% Momentum dependent width
if (m_daughter > 0)
    GAMMA = GAMMA0 * (m0./m) .* ((m.^2 - 4*m_daughter^2)./(m0^2 - 4*m_daughter^2)).^(3/2);
else
    GAMMA = GAMMA0;
end

f = 1 ./ (m - m0 + 1i*0.5*GAMMA); % Non-relativistic

end
