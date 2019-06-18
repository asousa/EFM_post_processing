close all; clear all;

% Load Campbell data
Mplate = 89.79;
% Csite = .105;
%tableVolts = [4.5,4,2, 1,.5,.2,0.1,0,-0.1, -.2,-.5,-1,-2,-4,-4.5];
% campbell_table_volts = [4482,3947,2033,987,504,208,102,0,-112,-203,-472,-994,-2005,-3897,-4545];
% efm_table_volts = [4620, 3985, 2097, 1043, 506, 187,  98, 0, -84, -208, -494, -1004, -2054, -3942, -4584];
% Load calibration voltages and vrefs:
voltages; 
ADC_SAMPLING_FREQ = 1000;
EFM_name = 'EFM011';
offset = 2;     % Samples to shift phase by (some mills have the optical encoder mis-aligned)
up_factor=10;   % Upsampling factor for peak detection (smooths out noise due to undersampling of sine wave)

output_filepath = "/Volumes/lairdata/EFM/Field Mill Post-Campaign Calibration/EFM calibration maps 5-12-2019";
input_filepath = "/Volumes/lairdata/EFM/Field Mill Post-Campaign Calibration/Calibration Data 5-12-2019";
prev_data_filepath = "/Volumes/lairdata/EFM/Field Mill Post-Campaign Calibration/Original calibration data/";
prev_campbell_filename = "Campbell Calib 080518-12am Map.mat";
% campbell_map_filename = "CR1000_EFM_arbitrary_voltages Map.mat";
campbell_map_filename = "CR1000_EFM_redo Map.mat";

gen_new_campbell_map = false;

table_separation = 9.25*0.054; % meters
%% Setup the Import Options
% if(strcmp(questdlg('Generate New Campbell Map?','Campbell Calib'),'Yes'))
if(gen_new_campbell_map)
    [filename, pathname] = uigetfile('.dat','Campbell raw data file');
    opts = delimitedTextImportOptions("NumVariables", 3);

    % Specify range and delimiter
    opts.DataLines = [5, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["TIMESTAMP", "RECORD", "E_field"];
    opts.VariableTypes = ["string", "double", "double"];
    opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
    opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Import the data
    tbl = readtable(filename, opts);
    TIMESTAMP = tbl.TIMESTAMP;
    startinds = find(tbl.RECORD==1);
    RECORD = tbl.RECORD(startinds(end):end);
    E_fieldC = tbl.E_field(startinds(end):end);

    campbell_field = zeros(1,length(campbell_table_volts));
    h=figure('Name','Campbell Calib File');
    set(gcf, 'Position', [0 100 595*2.5 421*2]);
    plot(E_fieldC)
    for ii=1:length(campbell_table_volts)
        waitfor(msgbox(['Please select on the graph the points corresponding to the ' ...
            'beginning and end of the E_field step corresponding to ' ...
            num2str(campbell_table_volts(ii)) ' kV/m'],'Calibrati`on Step Selection'))
        [xind,~] = ginput(2);
        campbell_field(ii) = mean(E_fieldC(round(xind(1):xind(2))))*Mplate;
    end
    close(h)
    
    save([filename(1:end-4) ' Map.mat'],'campbell_table_volts','campbell_field')

else
    % Load campbell calib map file
%     [filename, pathname] = uigetfile('.mat','Campbell Calib Map File');
%     load(fullfile(pathname,filename));
    load(fullfile(input_filepath, campbell_map_filename));
end


%% Load EFM File
efm_table_volts = EFM_table_volts(EFM_name);
ADC_REF = EFM_Vref(EFM_name);

filename = fullfile(input_filepath, [EFM_name '.bin']);
fileID = fopen(filename, 'r');
% [filename, pathname] = uigetfile('.bin');
% fileID = fopen(fullfile(pathname,filename), 'r');
data = fread(fileID,[1,3600*ADC_SAMPLING_FREQ],'uint16','n'); %% Change 2 to 1 in newer version
fclose(fileID);

data(2,:) = bitget(data(1,:),ones(1,length(data))*1,'uint16'); % Newer version injects bit for phase
data(1,:) = bitset(data(1,:),ones(1,length(data))*1,ones(1,length(data))*0,'uint16');
data = transpose(data);
time = 1/(ADC_SAMPLING_FREQ)*(0:3600*(ADC_SAMPLING_FREQ)-1);
time = transpose(time);

% Extract Amplitude Data

time = time(data(:,1)~=0);
sig = data((data(:,1)~=0),1)/65535*ADC_REF;
sig = detrend(sig,'linear',1:1e7:length(sig));

sig = lowpass(sig,110,ADC_SAMPLING_FREQ);

[su, tu] = resample(sig, time, up_factor*ADC_SAMPLING_FREQ);

[~, pklocpos] = findpeaks(su);
[~, pklocneg] = findpeaks(-su);
[envu,envl] = envPeak(su,pklocpos,pklocneg);
mag_up = (envu-envl)/2;
[mag, time] = resample(mag_up, tu, ADC_SAMPLING_FREQ);

% % Find pks will be used in phase as well. Run once here for efficiency
% [~, pklocpos] = findpeaks(sig);
% [~, pklocneg] = findpeaks(-sig);
% 
% [envu,envl] = envPeak(sig,pklocpos,pklocneg);
% % sighb = hilbert(sig);
% % mag = abs(sighb);
% mag = (envu-envl)/2;
E_field = mag; % No polarity yet

%% Extract Phase Data

phase = data((data(:,1)~=0),2);

phase = 2.*phase - 1;  % normalize to plus/minus 1

% phase = lowpass(phase,50,ADC_SAMPLING_FREQ);
% phase = sign(phase);
% Roll the phase left or right, if needed
phase = circshift(phase,offset); 

% upsample the phase vector as well
[pu, tu] = resample(phase, time, up_factor*ADC_SAMPLING_FREQ);


% t_pkloc = [time(pklocpos);time(pklocneg)];
% phase_pkloc = [(phase(pklocpos)>=0)*1+(phase(pklocpos)<0)*(-1); ...
%     (-phase(pklocneg)>=0)*1+(-phase(pklocneg)<0)*(-1)];
t_pkloc = [tu(pklocpos);tu(pklocneg)];
phase_pkloc = [(pu(pklocpos)>=0)*1+(pu(pklocpos)<0)*(-1); ...
    (-pu(pklocneg)>=0)*1+(-pu(pklocneg)<0)*(-1)];
pol = interp1(t_pkloc,phase_pkloc,time,'linear','extrap');
pol = sign(pol);
% pol = ones(size(pol));

% Apply Polarity and Calibration
E_field = E_field.*(pol);
% 
a = 120000; b = 120100;
figure(1);

plot(sig(a:b)); hold on; plot(phase(a:b)); plot(pol(a:b))
legend('signal','phase','polarity')
ylim([-1.2,1.2])
%% Get Calibration Values from Plot
EFM_magnitudes = zeros(1,length(efm_table_volts));
%efmVolts = zeros(1,length(tableVolts));
h=figure('Name','EFM Calib File');
set(gcf, 'Position', [0 100 595*2.5 421*2]);
plot(E_field)
for ii=1:length(efm_table_volts)
    waitfor(msgbox(['Please select on the graph the points corresponding to the ' ...
        'beginning and end of the E_field step corresponding to ' ...
        num2str(efm_table_volts(ii)) ' kV/m'],'Calibration Step Selection'))
    [xind,~] = ginput(2);
    EFM_magnitudes(ii) = mean(E_field(round(xind(1):xind(2))));
end
close(h)
%% Load previous Campbell data
prev_campbell_data = load(fullfile(prev_data_filepath, prev_campbell_filename));

% Load previous calibration data
% prev_data = load(fullfile(prev_data_filepath, [EFM_name ' Map.mat']));
%% Interpolate campbell and EFM data, and generate the calibration curve

% Interpolate all onto the same table volts:
interp_table_volts = linspace(-5000, 5000, 100);
campbell_interp = interp1(campbell_table_volts, campbell_field,   interp_table_volts,'linear','extrap');
EFM_volts_interp = interp1(efm_table_volts, EFM_magnitudes, interp_table_volts,'linear','extrap');

% Interpolate the campbell data onto uniform EFM volt axis
efmVolts = linspace(-1, 1, 100);
E_field_calib  = interp1(EFM_volts_interp, campbell_interp, efmVolts,'linear','extrap'); 

% Save the calibrated map:
outfile = fullfile(output_filepath, EFM_name + "_map_"+ datestr(now,'YYYY-mm-dd'));
save(outfile,'efmVolts','E_field_calib')

%% Interpolate previous data onto the same axis, for plotting
prev_camp_interp = interp1(1000*prev_campbell_data.tableVolts,  prev_campbell_data.E_field_True, interp_table_volts,'linear','extrap');
ideal_field = -interp_table_volts/table_separation; % volts per meter!

prev_E_field   = prev_data.E_field_True;
prev_EFM_volts = (prev_data.efmVolts).*sign(prev_E_field); % Prev cal doesn't include sign
prev_E_field_calib = interp1(prev_EFM_volts, prev_E_field, efmVolts,'linear','extrap');

% %% Plot new vs old Campbell map:
% figure(1);
% plot(interp_table_volts, prev_camp_interp,interp_table_volts, campbell_interp, interp_table_volts,ideal_field);
% legend('Previous Campbell map', 'New Campbell map','Ideal field')
% xlabel('Table voltage');
% ylabel('Electric field');

%% Let's look at the drift in gain!

P = polyfit(efmVolts,E_field_calib,1);
new_gain = P(1);
new_offset = P(2);
new_fit = P(1)*efmVolts + P(2);


P = polyfit(efmVolts,prev_E_field_calib,1);
old_gain = P(1);
old_offset = P(2);
old_fit = P(1)*efmVolts + P(2);

% figure(3);
% plot(efmVolts, prev_E_field_calib, efmVolts, old_fit);

gain_error = 100*(new_gain - old_gain)/new_gain;
offset_error = 100*(new_offset - old_offset)/new_offset;
fprintf('\nGain error:\t%6.2f percent\nOffset error:\t%6.2f percent\n',gain_error, offset_error);

%% Plot a nice figure, with new and old calibration curves, and linear fits:
f = figure();
set(groot,'defaultfigurecolor',[1 1 1])
set(groot,'defaultAxesFontSize',16)
set(groot,'defaultTextFontSize',18)
set(groot,'defaultAxesFontWeight','bold');
set(groot,'defaultTextFontWeight','bold');
set(groot,'defaultAxesLineWidth',2);
set(groot,'defaultUicontrolFontName','Arial');
set(groot,'defaultUitableFontName','Arial');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultTextFontName','Arial');
set(groot,'defaultUipanelFontName','Arial');
set(gcf, 'Position', [0 100 595*1 425*1.0]);

% hold on;
% p = plot(efmVolts, E_field_calib/1000,efmVolts, prev_E_field_calib/1000, 'LineWidth',2);
p = plot(efmVolts, E_field_calib/1000, 'LineWidth',2);

xlim([-1,1]);
ylim([-25,25]);
leg = legend('May 2019');%,'August 2018');
xlabel('EFM value (V)');
title([EFM_name ' ' datestr(now,'YYYY-mm-dd')]);
ylabel('Electric field (kV/m)');
set(leg, 'Location', 'northwest');
grid on;

new_str = sprintf('E_{new} = %3.3f V_{EFM} + %3.3f (kV/m)',new_gain/1000, new_offset/1000);
% prev_str = sprintf('E_{prev} = %3.3f V_{EFM} + %3.3f (kV/m)',old_gain/1000, old_offset/1000);
% drift_str = sprintf('Gain change:      %2.2f %%\nOffset change:   %2.2f %%', gain_error, offset_error);
text(0,-15, new_str,'FontSize',13);
% text(0,-18, prev_str,'FontSize',13);
% text(0,-22, drift_str,'FontSize',13);
outfile = fullfile(output_filepath, EFM_name + "_map_"+ datestr(now,'YYYY-mm-dd') + '.png');
saveas(f, outfile);