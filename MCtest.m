% Test sampling of Breit-Wigner amplitudes
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;

M0 = 1.5;
Gamma0 = 0.25;

% As a function of m, where sqrt(2m) from the jacobian at absolute squared
% level
g = @(m)  sqrt(2*m) ./ (m^2 - M0^2 + 1i*M0*Gamma0);

% As a function of m^2
G = @(m2) 1 ./ (m2 - M0^2 + 1i*M0*Gamma0);

N = 1e6;

gMAX = abs(g(M0))^2;
GMAX = abs(G(M0^2))^2;

g_sample = [];
G_sample = [];
a = 1.0; b = 2;

for i = 1:N
    m  = a + (b-a) * rand(1);       % flat sampling m
    g_ = abs(g(m))^2;
    
    if (g_ > rand(1)*gMAX)
        g_sample(end+1) = m;
    end
end
for i = 1:N
    m2 = a^2 + (b^2 - a^2) * rand(1); % flat sampling in m^2
    G_ = abs(G(m2))^2;
    
    if (G_ > rand(1)*GMAX)
        G_sample(end+1) = sqrt(m2);
    end
end


%%
close all;

[a,b] = hist(g_sample,200);
plot(b,a / sum(a), 'r-'); hold on;
[a,b] = hist(G_sample,200);
plot(b,a / sum(a), 'b-'); hold on;

axis tight;
xlabel('$M$ (GeV)','interpreter','latex');
ylabel('normalized counts','interpreter','latex');
l = legend('$g(M)$','$G(M^2)$'); set(l,'interpreter','latex');
axis([1 2 0 1.25*max(a/sum(a))]);
