% Total fit function
%
% mikael.mieskolainen@cern.ch, 2019

function [a2, y_tot, y, bg] = totfunc(XVAL, param, K, BWmode)

global m_daughter;

XVAL = XVAL(:)';
[A,m0,GAMMA0,phi, A_BG,gamma_BG,D_BG] = vec2param(param, K);

% Evaluate the BW function
y = cell(K,1);
for k = 1:K
    
    % Choose which Breit-Wigner parametrization we use
    mass = 0;
    
    if (BWmode(k) == 4 || BWmode(k) == 44)
        if (BWmode == 44)
            mass = m_daughter;
        end
        y{k} = ( A(k)*exp(1i*phi(k)) ) *    relativisticBWamplitude(XVAL, m0(k), GAMMA0(k), mass);
    end
    if (BWmode(k) == 5 || BWmode(k) == 55)
        if (BWmode == 55)
            mass = m_daughter;
        end
        y{k} = ( A(k)*exp(1i*phi(k)) ) * nonrelativisticBWamplitude(XVAL, m0(k), GAMMA0(k), mass);
    end
end

% Continuum (background)
bg = 1i * A_BG * BG1(XVAL, gamma_BG, D_BG);

% Total complex function
y_tot = zeros(1,length(XVAL));
for k = 1:length(y)
    y_tot = y_tot + y{k};
end
y_tot = y_tot + bg;

% Amplitude squared
a2 = abs(y_tot).^2 ; 

end