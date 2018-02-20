% Matlab m-file for ECE 445 NMOS Analysis

% If the HspiceToolbox is not in your path, uncomment and modify the following line
%addpath('/users/Kotecki/CppSim/CppSimShared/HspiceToolbox');
clear variables;
hspc_filename = 'NMOS.hspc';
% generate_hspc(hspc_filename, 'nmos_l', 'nmos_w', 'nmos_m', 'vgs', 'vds');

%% Set paramaters: vds, nmos_l, nmos_w and m
vds = 1.2;
nmos_l = 65e-9;
nmos_w = 100e-9;
nmos_m = 1;

%% Write paramaters and NGspice control statement
hspc_set_param('vds', vds, hspc_filename);
hspc_set_param('nmos_l', nmos_l, hspc_filename);
hspc_set_param('nmos_w', nmos_w, hspc_filename);
hspc_set_param('nmos_m', nmos_m, hspc_filename);

hspc_addline('.dc Vds 0 ''vds'' .01', hspc_filename); % DC sweep from 0 to vds

%% Figure size and location
Fig1 = figure('Position', [200, 75, 850, 600]);
grid on;
hold on;

%% Loop and run NGspice to generate IDS vs VDS for values of VGS
for vgs = .2:.2:vds
    hspc_set_param('vgs', vgs, hspc_filename); % set value of vgs
    legendname = sprintf('Vgs = %0.1fV',vgs); % define legend name
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

%% end of M file
