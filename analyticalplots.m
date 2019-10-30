% Analytical lineshape / Breit-Wigner plots for rho(770) case
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;
addpath ./src/

% Rho(770)
m0    = 0.775; % GeV
GAMMA = 0.15;  % GeV
m     = linspace(m0 - 2*GAMMA, m0 + 2*GAMMA, 1e4);

f1 = figure;
m_daughter = 0.139; % GeV

y = abs( nonrelativisticBWamplitude(m, m0,GAMMA, 0) ).^2;
plot(m,y/max(y),'b'); hold on;

% momentum dependent width
y = abs( nonrelativisticBWamplitude(m, m0,GAMMA, m_daughter) ).^2;
plot(m,y/max(y),'b--');

y = abs( relativisticBWamplitude(m, m0,GAMMA, 0) ).^2;
plot(m,y/max(y),'r');

% momentum dependent width
y = abs( relativisticBWamplitude(m, m0,GAMMA, m_daughter) ).^2;
plot(m,y/max(y),'r--');

plot(linspace(min(m),max(m),10), 0.5*ones(10,1), 'k:');
plot(m0*ones(1,10), linspace(0,1.1,10), 'k--');
plot((m0-GAMMA/2)*ones(1,10), linspace(0,1.1,10), 'k:');
plot((m0+GAMMA/2)*ones(1,10), linspace(0,1.1,10), 'k:');

axis([0.5 1.0 0 1.5]);

xlabel('$M$ (GeV)','interpreter','latex');
ylabel('$f(M)$ (a.u.)','interpreter','latex');
grid off;

l = legend('NR-BW', 'NR-BW (mom)', 'Rel. BW', 'Rel. BW (mom)'); legend('boxoff');

set(l,'interpreter','latex','fontsize',8);
axis tight; axis square;

filename = sprintf('../figs/breitwigners.pdf');
print(f1, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));





