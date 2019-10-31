% Fitting Central Exclusive Production Spectrum
%
% mikael.mieskolainen@cern.ch, 2019

% We need some histogramming functions (see other packages)
addpath /home/user/cernbox/#matlabcodes/

% Local functions
addpath ./src

clear; close all;


global XVAL;       % DATA
global YVAL;
global MDATA;      % Events
global FITRANGE;   % Fit range (min,max)
global K;          % Number of breit-wigners
global BWmode;     % Choose particular BW function
global massthresh; % Two particle mass threshold
global m_daughter; % Decay daughter mass

global CONSTRAINTS; % Enforce constraints
CONSTRAINTS = true;

modelist = 1:2;


for mode = modelist
    
    if (mode == 1)
        
        % x-ticks
        xticking = 0.2:0.2:2.6;
        
        particlemode = 'pipi';
        massthresh = 2*0.13957018; % 2xPion mass
        m_daughter = massthresh/2;
        
        latexmode = '\pi^{+}\pi^{-}';
        MDATA = csvread('masses_pipi.txt');
        
        FITRANGE = [massthresh, inf]; % GeV
        
        [YVAL,XVAL] = hist(MDATA, 220);
        
        BWmode = [44 44 44 44];
        K = length(BWmode);
        
        A      = [0.5 0.5 2.0 0.3];
        m0     = [0.77 0.97 1.275 1.52]; % GeV
        GAMMA0 = [0.15 0.08 0.185 0.1]; % GeV
        
        phi      = [pi 4/3*pi pi 2*pi];    % Radians
        A_BG     = 40;
        gamma_BG = 0.27;
        D_BG     = 0.19;
    end
    
    if (mode == 2)
        
        % x-ticks
        xticking = 0.95:0.2:2.55;
        
        particlemode = 'KK';
        massthresh = 2*0.493677; % 2xKaon mass
        m_daughter = massthresh/2;
        
        latexmode = 'K^{+}K^{-}';
        MDATA = csvread('masses_KK.txt');
        
        FITRANGE = [1.01, inf]; % GeV
        
        [YVAL,XVAL] = hist(MDATA, 80);
        
        BWmode = [44 44];
        K = length(BWmode);
        
        A      = [0.08 0.05];
        m0     = [1.542 1.732];    % GeV
        GAMMA0 = [0.078 0.123];   % GeV
        
        phi    = [pi pi/2];    % Radians
        A_BG   = 15.0;
        gamma_BG = 0.27;
        D_BG   = 0.4;
    end
    
    % PDG reference values
    global PDG_m0;
    global PDG_GAMMA0;
    global PDG_m0_err;
    global PDG_GAMMA0_err;
    
    if (mode == 1)
        PDG_name       = {'$\rho(770)$', '$f_0(980)$', '$f_2(1270)$', '$f_0(1500)$'};
        PDG_m0         = [0.77526 0.991 1.2751 1.505];
        PDG_GAMMA0     = [0.1491  0.075 0.1851 0.109]; 
        PDG_m0_err     = [0.25e-3 20e-3 1.2e-3 6e-3];
        PDG_GAMMA0_err = [0.8e-3  25e-3 2.9e-3 7e-3];
    end
    
    if (mode == 2)
        PDG_name       = {'$f_2''(1525)$', '$f_0(1710)$'};
        PDG_m0         = [1.525 1.722];
        PDG_GAMMA0     = [0.073  0.135]; 
        PDG_m0_err     = [5e-3 6e-3];
        PDG_GAMMA0_err = [6e-3  7e-3]; 
    end


%% Fit data

urand = @(a,b,n) a + (b-a) * rand(n,1);

%%if (mode == 1)
best_x0 = [];
best_chi2 = 1e9;
    
x0 = param2vec(A,m0,GAMMA0,phi, A_BG,gamma_BG,D_BG);

for j = 1:30
    x0 = fminsearch(@costfunc1, x0, optimset('TolX',1e-8,'MaxFunEvals',1e3,'MaxIter',1e4) );

    [chi2_fmin,~,ndf] = costfunc1(x0);
    fprintf('fminsearch: Chi2/ndf = %0.9f / %0.0f = %0.2f \n', chi2_fmin, ndf, chi2_fmin/ndf);

    if (chi2_fmin < best_chi2)
        best_chi2 = chi2_fmin;
        best_x0 = x0;
    end
end

%% Pick the best solution
x0 = best_x0;

% The best solution with constraint cost
CONSTRAINTS = false;
[chi2_fmin,~,ndf] = costfunc1(x0);
CONSTRAINTS = true;


%% Refine the solution and get errors

% Covariance matrix is: CovB = pinv(J'*J)*MSE
% NOTE that the J already includes the data error covariance matrix C
% by the documentation of nlinfit, the function absorbs this to J
%
% inv(J'C*J) -> inv(J'*J)
%
% Poisson errors: sigma = sqrt(N) <-> sigma^2 = N
sigma2 = YVAL;
W = 1./sigma2; W(isinf(W)) = 1e-9;
options = statset('maxiter', 1000);

[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(XVAL,YVAL,@nlinfitwrapper,x0,options,'Weights',W);

[chi2_nlinf,~,ndf] = costfunc1(beta);
fprintf('nlinfit: Chi2/ndf = %0.9f / %0.0f = %0.2f \n', chi2_nlinf, ndf, chi2_nlinf/ndf);

if (chi2_nlinf < chi2_fmin)
   x0 = beta; 
   chi2 = chi2_nlinf;
else
   chi2 = chi2_fmin; 
end

% CI function
alpha = 0.34; % 1 sigma
bCIw = nlparci(x0,R,J,alpha);

% Take the 1 sigma errors
x0_err = sqrt(diag(CovB));


%% Maximum likelihood fit

% for i = 1:10 % Recursively
%     x0 = fminsearch(@maxlikelihood, x0, optimset('TolX',1e-8,'MaxFunEvals',1e3,'MaxIter',1e4) );
%     [chi2,~,ndf] = costfunc1(x0);
%     fprintf('Chi2/ndf = %0.9f / %0.0f = %0.2f \n', chi2, ndf, chi2/ndf);
% end
% 
% m = linspace(0.28,3,1e4); plot(m, likelihood(m, x0, K));


%% Calculate the amplitude

close all;

% Mass values and axis
m = linspace(FITRANGE(1), 2.6, 5e4);

% Evaluate the fit function at mass value in m
[~, y_tot, y, bg] = totfunc(m, x0, K, BWmode);

% Get the parameters
[A,m0,GAMMA0,phi, A_BG,gamma_BG,D_BG] = vec2param(x0, K);

% Parameter errors
[A_err, m0_err, GAMMA0_err, phi_err, A_BG_err, gamma_BG_err, D_BG_err] = vec2param(x0_err, K);


%% PDG values and fit comparison

f1 = figure;

for k = 1:length(PDG_m0)
    %errbar(PDG_m0(k), PDG_GAMMA0(k), PDG_m0_err(k), 'k-', 'horiz');
    %errbar(PDG_m0(k), PDG_GAMMA0(k), PDG_GAMMA0_err(k), 'k-');
    
    rectangle('Position',[PDG_m0(k)-PDG_m0_err(k) PDG_GAMMA0(k)-PDG_GAMMA0_err(k) 2*PDG_m0_err(k) 2*PDG_GAMMA0_err(k)], ...
        'FaceColor',[0.9 0.9 0.9]); hold on;
    plot(PDG_m0(k), PDG_GAMMA0(k), 'k.', 'markersize', 7); hold on;
    
    plot(m0(k), GAMMA0(k), 'rs','markersize',2); hold on;
    rectangle('Position',[m0(k)-m0_err(k) GAMMA0(k)-GAMMA0_err(k) 2*m0_err(k) 2*GAMMA0_err(k)], ...
        'EdgeColor',[0.9 0 0]);
    
    text(PDG_m0(k)+0.033, PDG_GAMMA0(k)-0.005, PDG_name{k}, 'interpreter','latex');
    
    %errbar(m0(k), GAMMA0(k), GAMMA0_err(k), 'r-');
    %errbar(m0(k), GAMMA0(k), m0_err(k), 'r-', 'horiz');
end

% Draw geometry
plot(PDG_m0, PDG_GAMMA0, 'k--');
plot(m0, GAMMA0,'r:');

% Draw cross lines
if (mode == 1)
    vv = [1 3; 2 4];
    for k = 1:2
        plot(PDG_m0(vv(k,:)), PDG_GAMMA0(vv(k,:)), 'k--');
        plot(m0(vv(k,:)), GAMMA0(vv(k,:)), 'r:');
    end
end

xlabel('Mass $M_0$ (GeV)','interpreter','latex');
ylabel('Width $\Gamma_0$ (GeV)', 'interpreter','latex');
l = legend('PDG','Data'); set(l,'interpreter','latex'); legend('boxoff');
axis square;

if (mode == 1)
    axis([0.6 1.8 0.04 0.22]);
end
if (mode == 2)
    axis([1.5 1.85 0.0 0.3]);
end

filename = sprintf('../figs/%s_PDG.pdf', particlemode);
print(f1, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));


%% Mass spectrum together with data

f2 = figure;

ff = abs(y_tot).^2;

% Prediction intervals (NOT TESTED)
% funcwrap = @(beta,x) totfunc(x,beta,K);
% [Ypred,delta] = nlpredci(funcwrap, m, beta, R, 'Covar', CovB);
% plotfill(m, ff + delta, ff - delta, [0 0 1], [0 0 1], 0.05); hold on;

% Total
plot(m, ff, 'k-', 'linewidth', 0.8); hold on;

% How much space on top
if (mode == 1)
   yaxismax = prctile(ff, 99.99)*1.5;
elseif (mode == 2)
   yaxismax = prctile(ff, 99.99)*1.7;
end

cc = [0.5 0 0];
plot(XVAL,YVAL,'.','color',cc,'markersize',7);
errorbar(XVAL,YVAL,sqrt(YVAL), '.', 'color', cc, 'CapSize',0);

%title(sprintf('$\\chi^2$/ndf = %0.1f/%0.0f = %0.2f', chi2, ndf, chi2/ndf),'interpreter','latex');
xlabel(sprintf('$M_{%s}$ (GeV)', latexmode),'interpreter','latex');
ylabel(sprintf('Events / (%0.3f GeV)', XVAL(2)-XVAL(1)),'interpreter','latex');
axis square;
axis([min(xticking) max(m) 0 yaxismax]);
set(gca,'XTick',xticking)
set(gca,'yscale','log');

l = legend('$|f|^2$','Data');
set(l,'interpreter','latex'); legend('boxoff');

%errorbarwidth(h,100); % Must be last command

filename = sprintf('../figs/%s_complexfit0.pdf', particlemode);
print(f2, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));


%% Mass spectrum

f3 = figure;

% Total
plot(m, abs(y_tot).^2, 'k-', 'linewidth', 0.8); hold on;

% Individual resonances
for k = 1:length(y)
    plot(m, abs(y{k}).^2);
end

% Background/continuum
plot(m, abs(bg).^2, 'r--');

% Without background
plot(m, abs(y_tot - bg).^2);

xlabel(sprintf('$M_{%s}$ (GeV)', latexmode),'interpreter','latex');
ylabel(sprintf('Events / (%0.3f GeV)', XVAL(2)-XVAL(1)),'interpreter','latex');

axis([min(xticking) max(m) 0 yaxismax]);

axis square;
set(gca,'XTick',xticking)

% Use mod(radians, 2pi) to warp the phases
legcell = cell(length(y) + 4, 1);
for i = 1:length(y)
    legcell{i+1} = sprintf('$|r_{%d}|^2 \\; (M_0,\\Gamma_0)=(%0.3f\\pm%0.3f,%0.3f\\pm%0.3f), \\, (A,\\phi)=(%0.1f\\pm%0.1f, \\,%0.0f\\pm%0.0f)$', ...
        i, m0(i), m0_err(i), GAMMA0(i), GAMMA0_err(i), A(i), A_err(i), rad2deg(mod(phi(i),2*pi)), rad2deg(phi_err(i)));
end
legcell{1}     = '$|f|^2$';
legcell{end-1} = '$|f - \psi|^2$';
legcell{end-2} = sprintf('$|\\psi|^2 \\; A = %0.1f\\pm%0.1f, \\gamma = %0.2f\\pm%0.2f, D = %0.2f\\pm%0.2f$', ...
A_BG, A_BG_err, gamma_BG, gamma_BG_err, D_BG, D_BG_err);
legcell{end}   = sprintf('Data, $N = %d$', sum(YVAL));
title(sprintf('$\\chi^2$/ndf = %0.1f/%0.0f = %0.2f', chi2, ndf, chi2/ndf),'interpreter','latex');

%DATA
%h = errorbar(XVAL(1:jump:end),YVAL(1:jump:end),sqrt(YVAL(1:jump:end)), '.', 'color', [0.4 0 0]);
cc = [0.5 0 0];
plot(XVAL,YVAL,'.','color',cc,'markersize',7);
errorbar(XVAL,YVAL,sqrt(YVAL), '.', 'color', cc, 'CapSize',0);

l = legend(legcell);
set(l, 'interpreter','latex','fontsize',6);
legend('boxoff');

%errorbarwidth(h,80); % Must be last command

filename = sprintf('../figs/%s_complexfit.pdf', particlemode);
print(f3, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));


%% Difference and ratio

f4 = figure;
[~,YHAT] = costfunc1(x0);

subplot(1,2,1);
stepfill(XVAL, sqrt(YVAL), -sqrt(YVAL), [1 0 0], [1 0 0], 0.075); hold on;
plot(linspace(0,max(m), 10), zeros(10,1), 'k:');
stephist(XVAL, YVAL - YHAT, 'k-');

axis([min(xticking) max(m) -std(sqrt(YVAL))*6 std(sqrt(YVAL))*6]);
%set(gca,'XTick',xticking);
axis square;
ylabel('Data - Fit (Events)','interpreter','latex');
xlabel(sprintf('$M_{%s}$ (GeV)', latexmode),'interpreter','latex');

subplot(1,2,2);
% Ratio error
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulas
ratio_error = sqrt( 2*(sqrt(YVAL)./YVAL).^2 );
ratio_error(isnan(ratio_error)) = 0;
stepfill(XVAL, ratio_error+1, 1-ratio_error, [1 0 0], [1 0 0], 0.075);  hold on;
plot(linspace(0,max(m), 10), ones(10,1), 'k:'); hold on;
stephist(XVAL, YVAL./YHAT, 'k-');

axis([min(xticking) max(m) 0.25 1.75]);
%set(gca,'XTick',xticking);
axis square;
ylabel('Data / Fit','interpreter','latex');
xlabel(sprintf('$M_{%s}$ (GeV)', latexmode),'interpreter','latex');

filename = sprintf('../figs/%s_fitpull.pdf', particlemode);
print(f4, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));


%% Evaluate complex energy (mass) plane

mm = linspace(-1, 3, 1e3);
Z = zeros(length(mm));
Y = zeros(length(mm));

% Calculate values in the complex plane
for i = 1:length(mm)
    [a2, y_tot_, y_, bg_] = totfunc(mm(i) + 1i*mm, x0, K, BWmode);
    Z(i,:) = abs(y_tot_ - bg_).^2;
    Y(i,:) = abs(y_tot_).^2;
end

for i = 1:2
    fx = figure;
    if (i == 1)
        imagesc(mm, mm, log10(Z)'); hold on;
            title('log$_{10}$ $|f - b|^2$','interpreter','latex');
    elseif (i == 2)
        imagesc(mm, mm, log10(Y)'); hold on;
            title('log$_{10}$ $|f|^2$','interpreter','latex');
    end
    plot(linspace(min(mm),max(mm), 10), zeros(10,1), 'k:');
    plot(zeros(10,1), linspace(min(mm),max(mm), 10), 'k:');

    axis square; colormap(hot(500)); colorbar;
    xlabel(sprintf('Re[$M]_{%s}$ (GeV)', latexmode),'interpreter','latex');
    ylabel(sprintf('Im[$M]_{%s}$ (GeV)', latexmode),'interpreter','latex');
    
    set(gca,'YDir','normal')
    
    filename = sprintf('../figs/%s_complex_mass_%d.pdf', particlemode, i);
    print(fx, filename, '-dpdf');
    system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));
end


%% Derivative of the mass spectrum ~ "wavelet"

f5 = figure;
plot(m(1:end-1), diff(abs(y_tot).^2), 'r-'); hold on;
plot(linspace(0,max(m), 10), zeros(10,1), 'k:');
axis([min(xticking) max(m) -inf inf]);
axis square;
title('Derivative spectrum','interpreter','latex');
set(gca,'XTick',xticking);
l = legend('$\frac{\partial}{\partial M} |f(M)|^2$');
set(l, 'interpreter','latex'); legend('boxoff');

xlabel(sprintf('$M_{%s}$ (GeV)', latexmode),'interpreter','latex');
ylabel('(a.u.)','interpreter','latex');

filename = sprintf('../figs/%s_diff.pdf', particlemode);
print(f5, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));


%% Arc-length plots

% Integrate arg-length on the complex plane
[arclen, seglen] = arcintegral(real(y_tot), imag(y_tot));

f6 = figure;
subplot(1,2,1); % Arbitrary normalization by * sum(seglen)

lengths = seglen * sum(seglen);
[peakvals,peakind] = findpeaks(lengths);

% Remove spurious peak
if (mode == 1)
   peakvals = peakvals(1:end-1);
   peakind  = peakind(1:end-1);
end

plot(m(1:end-1), lengths, 'k-'); hold on;
legcell = cell(length(peakind)+1,1);
legcell{1} = '$l(M)$';

for k = 1:length(peakind)
   plot(m(peakind(k))*ones(10,1), linspace(0,peakvals(k),10));
   legcell{k+1} = sprintf('$r_%d = %0.3f$ GeV', k, m(peakind(k)));
end
xlabel(sprintf('$M_{%s}$ (GeV)', latexmode),'interpreter','latex');
ylabel('(a.u.)','interpreter','latex');
axis square;
title('$l(M) = \lim_{\Delta m\rightarrow 0} \int_M^{M+\Delta m} |f''(k)|\,dk$','interpreter','latex');
linspace(0,2.5);
axis([min(xticking) max(m) 0 inf]);
l = legend(legcell); set(l,'interpreter','latex','fontsize',7); legend('boxoff');

subplot(1,2,2);
derivative = diff(seglen * sum(seglen));
plot(m(1:end-2), derivative, 'r-'); hold on;
plot(linspace(0,max(m), 10), zeros(10,1), 'k:');
%plot(m(1:end-1), diff(abs(y_tot).^2), 'r-');
axis([min(xticking) max(m) -inf inf]);
axis square;
title('$\frac{\partial}{\partial M} l(M)$','interpreter','latex');

xlabel(sprintf('$M_{%s}$ (GeV)', latexmode),'interpreter','latex');
ylabel('(a.u.)','interpreter','latex');

filename = sprintf('../figs/%s_derivative.pdf', particlemode);
print(f6, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));


%% Trajectory in the C-plane

f7 = figure;
h = plot(real(y_tot), imag(y_tot), 'k-'); hold on;

% Create indices uniformly along segment lengths
Nind = 55;
steplength = arclen / Nind;
indices = [];
counter = 0;
for i = 1:length(seglen)
    counter = counter + seglen(i);
    if (counter >= steplength)
        indices(end+1) = i;
        counter = 0;
    end
end
% Plot mass values along the C-plane curve
plot(real(y_tot(indices)), imag(y_tot(indices)),'k.');

for i = indices
    text(real(y_tot(i)), imag(y_tot(i)), sprintf('%0.2f', m(i)), 'fontsize', 6);
end
% Scan the pole (resonance) indices
pole_index = zeros(length(y),1);
for k = 1:length(m0)
    [min_val, min_ind] = min(abs(m - m0(k)));
    pole_index(k) = min_ind;
end
% Plot poles on the C-plane
for k = 1:length(m0)
    i = pole_index(k);
    % Pole
    re_coord = real(y_tot(i));  
    im_coord = imag(y_tot(i));
    plot(re_coord, im_coord, 'ko', 'markersize', 7); hold on;
    
    % "Phasor"
    plot([0 re_coord], [0 im_coord], '-');
end

xlabel('Re $[f]$ (a.u.)','interpreter','latex');
ylabel('Im $[f]$ (a.u.)','interpreter','latex');
axis equal;

% Calculate the center of gravity
cg = [mean(real(y_tot(indices))), mean(imag(y_tot(indices)))];
cg2 = [std(real(y_tot(indices))), std(imag(y_tot(indices)))];

scale = 2.4;
xbound = [cg(1)-scale*cg2(1), cg(1)+scale*cg2(1)];
ybound = [cg(2)-scale*cg2(2), cg(2)+scale*cg2(2)];

% horizontal and vertical axis
plot(linspace(xbound(1), xbound(2), 3), zeros(3,1), 'k:');
plot(zeros(3,1), linspace(ybound(1), ybound(2),3), 'k:');

axis([xbound(1) xbound(2) ybound(1) ybound(2)]);
l = legend('$f(M) \in \mathbf{C}$',sprintf('$M_{%s}$ (GeV)', latexmode),'Resonance'); set(l,'interpreter','latex');
legend('boxoff');

filename = sprintf('../figs/%s_cplane.pdf', particlemode);
print(f7, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));


%% Create common legend

legcellc = cell(length(y) + 2, 1);
legcellc{1}     = '$f$';
for i = 1:length(y)
    legcellc{i+1} = sprintf('$r_{%d}$', i);
end
legcellc{end} = '$b$';


% Modulo-phase plots

% LEFT
f8 = figure;
%yyaxis left;
plot(m, abs(y_tot), 'k-'); hold on;
for k = 1:K
    plot(m, abs(y{k})); hold on;
end
plot(m, abs(bg),'r--'); hold on;
xlabel(sprintf('$M_{%s}$ (GeV)', latexmode),'interpreter','latex');
ylabel('$|f|$ (a.u.)','interpreter','latex');
title('Modulus','interpreter','latex');
axis square;
axis([min(xticking) max(m) 0 inf]);
set(gca,'XTick',xticking);
l = legend(legcellc); legend('boxoff'); set(l,'interpreter','latex');

filename = sprintf('../figs/%s_amplitude.pdf', particlemode);
print(f8, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));


% RIGHT
f9 = figure;
%yyaxis right;
%plot(m, angle(y_tot)); hold on;
%plot(m, atan2(imag(y_tot),real(y_tot))); hold on;
plot(m, rad2deg(unwrap(angle(y_tot))), 'k-'); hold on;

for k = 1:K
    plot(m, rad2deg(unwrap(angle(y{k}))) ); hold on;
end
plot(m, rad2deg(unwrap(angle(bg))), 'r--'); hold on;

xlabel(sprintf('$M_{%s}$ (GeV)', latexmode),'interpreter','latex');
ylabel('arg$(f)$ (deg)','interpreter','latex');
title('Phase (unwrapped)','interpreter','latex');
axis square;
axis([min(xticking) max(m) 0 inf]);
set(gca,'XTick',xticking);
l = legend(legcellc); legend('boxoff'); set(l,'interpreter','latex','location','southeast');

filename = sprintf('../figs/%s_phase.pdf', particlemode);
print(f9, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));


%% Real-Imag plot

f10 = figure;
plot(m, real(y_tot), 'k-'); hold on;
for k = 1:K
    plot(m, real(y{k}), '-');
end
plot(m, real(bg), 'r--');
plot(linspace(0,max(m),10), zeros(10,1), 'k:'); hold on; % x-axis
xlabel(sprintf('$M_{%s}$ (GeV)', latexmode),'interpreter','latex');
ylabel('(a.u.)','interpreter','latex');
l = legend('Re $[f]$','Im $[f]$'); legend('boxoff');
set(l,'interpreter','latex');
title('Real part','interpreter','latex');
axis square;
axis([min(xticking) max(m) -inf inf]);
l = legend(legcellc); legend('boxoff'); set(l,'interpreter','latex','location','southeast');
set(gca,'XTick',xticking);

filename = sprintf('../figs/%s_real.pdf', particlemode);
print(f10, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));


f11 = figure;
plot(m, imag(y_tot), 'k-'); hold on;
for k = 1:K
    plot(m, imag(y{k}), '-');
end
plot(m, imag(bg), 'r--');
plot(linspace(0,max(m),10), zeros(10,1), 'k:'); hold on; % x-axis
xlabel(sprintf('$M_{%s}$ (GeV)', latexmode),'interpreter','latex');
ylabel('(a.u.)','interpreter','latex');
l = legend('Re $[f]$','Im $[f]$'); legend('boxoff');
set(l,'interpreter','latex');
title('Imaginary part','interpreter','latex');
axis square;
axis([min(xticking) max(m) -inf inf]);
l = legend(legcellc); legend('boxoff'); set(l,'interpreter','latex');
set(gca,'XTick',xticking);

filename = sprintf('../figs/%s_imag.pdf', particlemode);
print(f11, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));


%% Pairwise interference
% |a+b|^2 = |a|^2 + |b|^2 + 2*Re[a*b]

f12 = figure;

legs = {};
for i = 1:K
    for j = i:K
        if (i ~= j)
            interference = 2*real(conj(y{i}).*y{j});
            plot(m, interference); hold on;
            legs{end+1} = sprintf('2Re[$r_%d^*r_%d$]', i, j);
        end
    end
end
plot(linspace(0,max(m),10), zeros(10,1), 'k:'); hold on; % x-axis

xlabel(sprintf('$M_{%s}$ (GeV)', latexmode),'interpreter','latex');
ylabel('Events','interpreter','latex');
l = legend(legs);
set(l,'interpreter','latex'); legend('boxoff');

% if (mode == 1)
%    set(l,'location','southeast');
% end

axis square;
title('Pairwise interference','interpreter','latex');
axis([min(xticking) max(m) -inf inf]);
set(gca,'XTick',xticking)

filename = sprintf('../figs/%s_pairwise_interference.pdf', particlemode);
print(f12, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));


f13 = figure;
legs = {};
for i = 1:K
    interference = 2*real(conj(y{i}).*bg);
    plot(m, interference); hold on;
    legs{end+1} = sprintf('2Re[$r_%d^*b$]', i);
end

plot(linspace(0,max(m),10), zeros(10,1), 'k:'); hold on; % x-axis
xlabel(sprintf('$M_{%s}$ (GeV)', latexmode),'interpreter','latex');
ylabel('Events','interpreter','latex');
l = legend(legs);
set(l,'interpreter','latex'); legend('boxoff');

% if (mode == 1)
%    set(l,'location','southeast');
% end

axis square;
title('Interference with continuum','interpreter','latex');
axis([min(xticking) max(m) -inf inf]);
set(gca,'XTick',xticking)

filename = sprintf('../figs/%s_continuum_interference.pdf', particlemode);
print(f13, filename, '-dpdf');
system(sprintf('pdfcrop --margins 2 %s %s', filename, filename));

end
