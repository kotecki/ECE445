% Matlab m-file for ECE 445 Inverter_Test2
% DEK (2/2018)

clear variables;
hspc_filename = 'Inverter_test2.hspc'; % Set hspc filename

%% Set paramaters:
vdd = 1.2;
mos_l = 65e-9;
nmos_w = 100e-9;
nmos_m = 1;
pmos_w = 200e-9;
pmos_m = 1;

%% Write paramaters
hspc_set_param('vdd', vdd, hspc_filename);
hspc_set_param('mos_l', mos_l, hspc_filename);
hspc_set_param('nmos_w', nmos_w, hspc_filename);
hspc_set_param('nmos_m', nmos_m, hspc_filename);
hspc_set_param('pmos_w', pmos_w, hspc_filename);
hspc_set_param('pmos_m', pmos_m, hspc_filename);

%% Set Control statement for DC Analysis an Run NGspice
hspc_addline('.dc V0 0 ''vdd'' .001', hspc_filename); % DC sweep from 0 to vdd
ngsim(hspc_filename);  % run ngspice

%% Load data and generate figure for DC response
data = loadsig('simrun.raw'); % load simulation results and extract vds and ids
Vin = evalsig(data,'vin'); % equally spaced input voltage
Vout1 = evalsig(data, 'vout1');
deriv = diff(Vout1)/(Vin(3) - Vin(2));

Fig1 = figure('Name', 'DC Characteristics', 'Position', [200, 75, 850, 600]);
lw = 2; % set linewidth
fs = 16; % set font size
plot(Vin, Vout1, 'linewidth', 1.5*lw) % plot Vout vs Vin
hold on;
V1 = linspace(0, vdd, 100); % generate Vout = Vin
plot(V1, V1, 'linewidth', lw); % plot Vout = Vin
grid on;
set(gca, 'fontsize', fs); % increase font size
xlabel('Vin (V)', 'fontsize', fs); % x-axis labels
ylabel('Vout (V)', 'fontsize', fs); % y-axis labels

%% Calculate intersection to determine the switching voltage Vs
[x0,y0] = intersections(Vin, Vout1, V1, V1, 0);
plot(x0, y0, 'k.', 'MarkerSize', 25); % plot Vs

%% Calculate and plot the noise margins
for index=1:length(deriv)-1 % estimagte locations of NML and NMH
    if (deriv(index+1) < -1 && deriv(index) > -1)
        NMLi = Vin(index);
        NMLo = Vout1(index);
    end
    if (deriv(index+1) > -1 && deriv(index) < -1)
        NMHi = Vin(index);
        NMHo = Vout1(index);
    end
end

color1 = 1/255*[162,20,70];
plot([NMLi, NMLi], [0, NMLo], '--', 'Color', color1, 'linewidth', lw); 
plot([NMLi-.1, NMLi+.1], [NMLo+.1, NMLo-.1], 'Color', color1, 'linewidth', lw);
plot([NMHi, NMHi], [0 NMHo], '--', 'Color', color1, 'linewidth', lw);
plot([NMHi-.1, NMHi+.1], [NMHo+.1, NMHo-.1], 'Color', color1, 'linewidth',lw);

if (NMLo + .08 > vdd || NMHo + .08 < 0)
    axis([-.1 vdd+.1 -.1 vdd+.1]);
else
    axis([0 vdd 0 vdd]);
end
     
%% Add Vs, NMH, and NML to graph
text1 = sprintf('V_S = %0.3g V \n', x0);
text2 = sprintf('NML = %0.3g V \n', NMLi);
text3 = sprintf('NMH = %0.3g V', vdd - NMHi);
text = [text1, text2, text3];
dim=[.6, .35, .2, .2];
annotation('textbox',dim,'String',text,'FitBoxToText','on', 'BackgroundColor','white', ...
    'FaceAlpha',0.5,'fontsize', fs);

%% Set Control statement for Transient Analysis an Run NGspice
hspc_addline('.tran 1e-12 5e-9 0 0', hspc_filename); % transient analysis
ngsim(hspc_filename);  % run ngspice

%% Load data and generate figure 
Fig2 = figure('Name', 'Transient Response', 'Position', [100, 75, 850, 600]);

data = loadsig('simrun.raw'); % load simulation results and extract vds and ids
time = evalsig(data,'TIME');
Vin = evalsig(data, 'vout3');
Vout = evalsig(data, 'vout4');
plot(1e9*time, Vin, 1e9*time, Vout, 'linewidth', lw) % plot Vout vs time
grid on;
legend('Vout3', 'Vout4');
set(gca, 'fontsize', fs); % increase font size
xlabel('Time (ns)', 'fontsize', fs); % x-axis labels
ylabel('Vout (V)', 'fontsize', fs); % y-axis labels
axis([0 1e9*max(time) -0.2 1.6]);

%% Calculate the propagation delay
vm = vdd / 2;

% tphl
Vin_index = [ ]; % Define empty array to calculate the vdd/2 crossings for Vin
Vout_index = [ ]; % Define empty array to calculate the vdd/2 crossings for Vout
for index = 1:length(Vin)-1 % calculate the vdd/2 crossing of Vin for the rising edge
    if(Vin(index) < vm && Vin(index+1) > vm)
        slope = (Vin(index+1) - Vin(index)) / (time(index + 1) - time(index));
        Vin_index(end+1) = time(index) - (Vin(index) - vm)./slope;
    end
end

for index = 1:length(Vout)-1 % calculate the vdd/2 crossings of Vout for the falling edge
    if(Vout(index) > vm && Vout(index+1) < vm)
        slope = (Vout(index+1) - Vout(index)) / (time(index + 1) - time(index));
        Vout_index(end+1) = time(index) - (Vout(index) - vm)./slope;
    end
end

if Vin_index(1) >= Vout_index(1)
    DeltaT = Vout_index(4) - Vin_index(3); % time shift
else
    DeltaT = Vout_index(4) - Vin_index(4); % time shift
end
tphl = DeltaT; % propagation delay low to high

% tplh
Vin_index = [ ]; % Define empty array to calculate the vdd/2 crossings for Vin
Vout_index = [ ]; % Define empty array to calculate the vdd/2 crossings for Vout
for index = 1:length(Vout)-1 % calculate the vdd/2 crossing of Vin for the rising edge
    if(Vout(index) < vm && Vout(index+1) > vm)
        slope = (Vout(index+1) - Vout(index)) / (time(index + 1) - time(index));
        Vout_index(end+1) = time(index) - (Vout(index) - vm)./slope;
    end
end

for index = 1:length(Vin)-1 % calculate the vdd/2 crossings of Vout for the falling edge
    if(Vin(index) > vm && Vin(index+1) < vm)
        slope = (Vin(index+1) - Vin(index)) / (time(index + 1) - time(index));
        Vin_index(end+1) = time(index) - (Vin(index) - vm)./slope;
    end
end

if Vin_index(1) >= Vout_index(1)
    DeltaT = Vout_index(4) - Vin_index(3); % time shift
else
    DeltaT = Vout_index(4) - Vin_index(4); % time shift
end
tplh = DeltaT; % propagation delay low to high

%% Add propagation delay to the graph
t0 = sprintf('Propagation Delay \n');
t1 = sprintf('t_{PLH} = %0.3g ps \n',tphl*1e12);
t2 = sprintf('t_{PHL} = %0.3g ps \n',tplh*1e12);
t3 = sprintf('t_P = %0.3g ps', (tphl + tplh)*1e12/2);
text=[t0, t1, t2, t3];
dim = [.15, .81 .1 .1];
annotation('textbox',dim,'String',text,'FitBoxToText','on', 'BackgroundColor','white', ...
    'FaceAlpha',0.5,'fontsize', fs);

%% Show DC characteristics on top
uistack(Fig1, 'top');

%% end of M file