% Matlab m-file for ECE 445 Ring Oscillator Test
% DEK (2/2018)

clear variables;
hspc_filename = 'RingOsc.hspc'; % Set hspc filename

%% Set paramaters:
nstages = 101; % number of stages in ring oscillator - used for calculation of t_p
vdd = 1.2;
mos_l = 65e-9;
nmos_w = 100e-9;
nmos_m = 1;
pmos_w = 200e-9;
pmos_m = 1;
temperature = 27; % temperature in C

%% Write paramaters
hspc_set_param('vdd', vdd, hspc_filename);
hspc_set_param('mos_l', mos_l, hspc_filename);
hspc_set_param('nmos_w', nmos_w, hspc_filename);
hspc_set_param('nmos_m', nmos_m, hspc_filename);
hspc_set_param('pmos_w', pmos_w, hspc_filename);
hspc_set_param('pmos_m', pmos_m, hspc_filename);

%% Set Control statement for Transient Analysis an Run NGspice
hspc_addline('.tran 1e-11 1e-8 0 2e-12', hspc_filename); % transient analysis
temp_string = sprintf('.temp %3g', temperature); % set temperature in hspc file
hspc_addline_continued(temp_string, hspc_filename); % set temperature
ngsim(hspc_filename);  % run ngspice

%% Load data and generate figure 
Fig1 = figure('Name', 'Transient Response', 'Position', [100, 75, 850, 600]);
lw = 2; % set linewidth
fs = 16; % set font size

data = loadsig('simrun.raw'); % load simulation results and extract vds and ids
time = evalsig(data,'TIME');
Vout = evalsig(data, 'vout');
plot(1e9*time, Vout, 'linewidth', lw) % plot Vout
grid on;
set(gca, 'fontsize', fs); % increase font size
xlabel('Time (ns)', 'fontsize', fs); % x-axis labels
ylabel('Vout (V)', 'fontsize', fs); % y-axis labels
title([num2str(nstages), ' Stage Ring Oscillator (',num2str(temperature),'^oC)']);

%% Calculate frequency using findpeaks
[pks, pks_i] = findpeaks(Vout, 'MinPeakDistance', 300); % locate peaks
pks_i(pks < vdd) = []; % elliminate peaks lower than vdd
pks(pks < vdd)= [];
hold on
plot(1e9*time(pks_i), pks, 'or');
legend('Vout', 'Overshoot');
fest = [];
for i = 1: length(pks)-1
    f_est(i) = 1 / (time(pks_i(i + 1)) - time(pks_i(i)));
end
f_mean = mean(f_est); % average frequency
f_std = std(f_est); % standard deviation of frequency
t_p = 1 / (f_mean * nstages * 2); 

%% Annotate graph
text0 = sprintf('Average Frequency (from peaks) = %0.4g MHz\n', 1e-6 * f_mean);
text1 = sprintf('Standard Deviation (from peaks) = %0.4g MHz\n', 1e-6 * f_std);
text2 = sprintf('Average Propagation Delay (1/2fn) = %0.4g ps', 1e12 * t_p);
text = [text0 text1 text2];
dim=[.15, .2, .2, .2];
annotation('textbox',dim,'String',text,'FitBoxToText','on', 'BackgroundColor','white', ...
    'FaceAlpha',0.85,'fontsize', fs);

%% end of M file