%% Matlab m-file for ECE 445 NMOS Analysis

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
Schematic_Name = sprintf('PMOS_65nm_test'); % name of schematic

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

fprintf(hspcfile, '**** Paramenter Statements ****\n');
vdd = -1.2;
fprintf(hspcfile, '.param pmos_l = 65e-9 \n');  
fprintf(hspcfile, '.param pmos_w = 100e-9 \n');
fprintf(hspcfile, '.param pmos_m = 10 \n');
fprintf(hspcfile, '.param vds = 0 \n');
fprintf(hspcfile, '.param vgs = 0\n\n');

fprintf(hspcfile, '**** Simulation Statement ****\n');
fprintf(hspcfile, '.dc Vds 0 ''vds'' -.01 \n\n');

fprintf(hspcfile, '**** Include Statements ****\n');
fprintf(hspcfile, '.include ../../../SpiceModels/ECE214_models.mod \n');
fprintf(hspcfile, '.include ../../../SpiceModels/ECE445_models.mod \n\n');

% fprintf(hspcfile, '**** Initial Conditions ****\n');
% fprintf(hspcfile, '.ic v(out1)=5 \n\n');

fprintf(hspcfile, '**** Simulation Options ****\n');
fprintf(hspcfile, '.options post=1 delmax=5p relv=1e-6 reli=1e-6 relmos=1e-6 method=gear \n');
fprintf(hspcfile, '.temp 27\n');
fprintf(hspcfile, '.global gnd \n');
fprintf(hspcfile, '.op \n\n');
fprintf(hspcfile, '**** End of NGspice hspc file \n');
fclose(hspcfile);

%% Section 3: Run simulation, Load results, and analyze data

Fig1 = figure('Position', [200, 75, 850, 600]);
grid on;
hold on;
hspc_set_param('vds', vdd, hspc_filename); % 

%% Loop and run NGspice to generate IDS vs VDS for values of VGS
for vgs = .0:-.2:vdd
    hspc_set_param('vgs', vgs, hspc_filename); % set value of vgs
    legendname = sprintf('Vgs = %0.2fV',vgs); % define legend name
    ngsim(hspc_filename);  % run ngspice
    data = loadsig('simrun.raw'); % load simulation results and extract vds and ids
    vds = evalsig(data,'VOLTAGE');
    ids = evalsig(data, 'i_vds');
    igs = evalsig(data, 'i_vgs');
    plot(vds, -1E3*ids, 'linewidth',2.0, 'displayname', legendname) % plot ids vs vds
end
hold off;
set(gca, 'fontsize', 16); % increase font size
legend('location', 'nw');
xlabel('Vds (V)', 'fontsize', 16); % x-axis labels
ylabel('Ids (mA)', 'fontsize', 16); % y-axis labels

%% end of .m file
