%% Matlab m-file for ECE 445 Ring Oscillator Test

%% Copyright (c) 2019 by David E. Kotecki. All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:

% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% Section 1: Define: CPPSim location, library, and schematic
clear variables;  
CppSim_Location = sprintf('C:/CppSim'); % location of CppSim directory
Design_Library = sprintf('ECE445_2019'); % name of design library
Schematic_Name = sprintf('RingOsc_test'); % name of schematic

%% Section 2: Generate HSPC file and run NGspice
addpath(sprintf('%s/CppSimShared/HspiceToolbox', CppSim_Location)); % add ngspice matlab toolbox to the path
Working_Dir = sprintf('%s/SimRuns/%s/%s', CppSim_Location, Design_Library, Schematic_Name);
if ~exist(Working_Dir, 'dir')  
    mkdir(Working_Dir)  % create working directory if it does not exist
end
cd(Working_Dir) % set current folder to the working directory
hspc_filename = sprintf('%s.hspc', Schematic_Name);   % define hspc filename

hspcfile = fopen(hspc_filename, 'w');    % open file 'hspc_filename' for write
fprintf(hspcfile, '**** NGspice HSPC file **** \n');
fprintf(hspcfile, '**** File: %s/%s **** \n', pwd, hspc_filename);
fprintf(hspcfile, '**** Date: %s **** \n\n', datestr(datetime('now')));

fprintf(hspcfile, '**** Simulation Statement ****\n');
fprintf(hspcfile, '.tran 5e-10 1.5e-8 0 6e-12 UIC \n\n');

fprintf(hspcfile, '**** Simulation Temperature ****\n');
temperature = 27;
fprintf(hspcfile, '.temp %3g \n', temperature);

fprintf(hspcfile, '**** Paramenter Statements ****\n');
nstages = 101; % number of stages - used for propagation delay calculation
vdd = 1.2;
fprintf(hspcfile, '.param vdd = %g \n', vdd);  
fprintf(hspcfile, '.param mos_l = 65e-9 \n');  
fprintf(hspcfile, '.param nmos_w = 100e-9 \n'); 
fprintf(hspcfile, '.param nmos_m = 1 \n');  
fprintf(hspcfile, '.param pmos_w = 200e-9 \n');  
fprintf(hspcfile, '.param pmos_m = 1 \n\n');

fprintf(hspcfile, '**** Include Statements ****\n');
fprintf(hspcfile, '.include ../../../SpiceModels/ECE214_models.mod \n');
fprintf(hspcfile, '.include ../../../SpiceModels/ECE445_models.mod \n\n');

fprintf(hspcfile, '**** Initial Conditions ****\n');
fprintf(hspcfile, '.ic v(vout)=0.0 \n\n');

fprintf(hspcfile, '**** Rail Voltage ****\n');
fprintf(hspcfile, 'vsup vdd 0 ''vdd''\n\n'); 

fprintf(hspcfile, '**** Simulation Options ****\n');
fprintf(hspcfile, '.options post=1 delmax=5p relv=1e-6 reli=1e-6 relmos=1e-6 method=gear \n');
fprintf(hspcfile, '.global gnd vdd \n');
fprintf(hspcfile, '.op \n\n');
fprintf(hspcfile, '**** End of NGspice hspc file ****\n');
fclose(hspcfile); 
 
ngsim(hspc_filename); % run ngspice  

%% Section 3: Load simulation results and analyze data
Fig1 = figure('Name', 'Transient Response', 'Position', [100, 75, 850, 600]);
lw = 2; % set linewidth
fs = 16; % set font size

data = loadsig('simrun.raw'); 
time = evalsig(data,'TIME');
Vout = evalsig(data, 'vout');
plot(1e9*time, Vout, 'linewidth', lw) % plot Vout
grid on;
set(gca, 'fontsize', fs); % increase font size
xlabel('Time (ns)', 'fontsize', fs); % x-axis labels
ylabel('Vout (V)', 'fontsize', fs); % y-axis labels
title([num2str(nstages), ' Stage Ring Oscillator (',num2str(temperature),'^oC)']);

%% Calculate frequency using findpeaks
[pks, pks_i] = findpeaks(Vout, 'MinPeakDistance', 200); % locate peaks
pks_i(pks < vdd) = []; % elliminate peaks lower than vdd
pks(pks < vdd)= [];
pks_i(1) = []; % elliminate first peak
pks(1) = [];
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

%% Plot power and calculate average power
Fig2 = figure('Name', 'Power', 'Position', [120, 75, 850, 600]);

Power = vdd .* evalsig(data, 'i_vsup');
skip = 100;
plot(1e9*time(skip:end), Power(skip:end), 'linewidth', lw) % plot power ignore first 50 points))
AveragePower = abs(mean(Power(skip:end)));
grid on;
set(gca, 'fontsize', fs); % increase font size
xlabel('Time (ns)', 'fontsize', fs); % x-axis labels
ylabel('Power (W)', 'fontsize', fs); % y-axis labels
title([num2str(nstages), ' Stage Ring Oscillator (',num2str(temperature),'^oC)']);
texta = sprintf('Average Power = %0.4g mW', AveragePower .*1e3);
dim=[.15, .35, .2, .2];
annotation('textbox',dim,'String',texta,'FitBoxToText','on', 'BackgroundColor','white', ...
    'FaceAlpha',0.85,'fontsize', fs);

%% Show transient characteristics on top
uistack(Fig1, 'top');
%% end of .m file
