% Compute local site corrections between Campbell and CU EFM measurements:
close all; clear all;
ADC_SAMPLING_FREQ = 1000;
raw_data_dir = "/Volumes/lairdata/EFM/RELAMPAGO Data/Campaign Data";
campbell_data_dir = "/Volumes/lairdata/EFM/RELAMPAGO Data/Campaign Data/Campbell Field Deployment";
cal_dir = "/Volumes/lairdata/EFM/Field Mill Post-Campaign Calibration/EFM calibration maps 6-17-2019";
fig_dir = "/Volumes/lairdata/EFM/Field Mill Post-Campaign Calibration/Site Corrections";

Mplate = 89.79;  % Campbell's mill-specific collector area correction.
Csite =  0.105;   % Campbell's reported site correction. Probably valid, but will vary with setup height

site_name = "Villa-del-Rosario"; EFM = 'EFM002'; lag = 0;  % works-ish
% site_name = "Villa-Carlos-Paz";  EFM = 'EFM008';  lag = 0;   % works-ish
% site_name = "Cordoba"; EFM = 'EFM011'; lag = 0;    % "CR1000_EFM_old.dat" kinda works... missing EFM data for the other file
% site_name = "Manfredi"; EFM = 'EFM004';  lag = 0;  % Raw data is garbage at calibration :/ 
% site_name = "Pilar";    EFM = 'EFM006';  lag = 10;  % This one works!

% site_name = "Almafuerte2"; EFM = 'EFM001'; lag = 30;

argentina_time_offset = hours(3); 
phase_offset = 1;
% Load Campbell data:
% campbell_file = fullfile(campbell_data_dir,EFM,'CR1000_EFM.dat');
% campbell_file = fullfile(campbell_data_dir,'RELAMPAGO Pilar Calib','Campbell-Pilar-121918','CR1000_EFM-cal-121918-1147am.dat');
campbell_file = fullfile(campbell_data_dir,'RELAMPAGO Pilar Calib','Campbell-Pilar-121918','CR1000_EFM-cal-121918-1153am.dat');
Cfile = readtable(campbell_file);
Cfile = Cfile(3:end,:);

timec = regexprep(Cfile.TIMESTAMP, ':\d\d$', '$&.0', 'lineanchors');
timec = datetime(timec, 'InputFormat', 'yyyy-MM-dd  HH:mm:ss.S');
timec = timec + argentina_time_offset;
EfieldC = str2double(Cfile.E_field);

EfieldC = EfieldC*Csite*Mplate;  % Apply site correction to Campbell
% Calibration map
cal_file = fullfile(cal_dir,sprintf("%s_map_2019-06-17.mat",EFM));  
cal_data = load(cal_file); % returns efmVolts, E_field_calib



%% Load corresponding EFM data:

% next_file_time = file_time + hours(1);
ndvec = datevec(timec(1) + hours(1));
next_file = fullfile(raw_data_dir,site_name,'DATA',...
            sprintf('%d',ndvec(1)),sprintf('%d',ndvec(2)), sprintf('%d',ndvec(3)),...
            sprintf('%02d.bin',ndvec(4)));
have_next   = isfile(next_file);

dvec = datevec(timec(1));
cur_file = fullfile(raw_data_dir,site_name,'DATA',...
            sprintf('%d',dvec(1)),sprintf('%d',dvec(2)), sprintf('%d',dvec(3)),...
            sprintf('%02d.bin',dvec(4)));

% Load current file
fprintf("Loading %s\n",cur_file);
fileID = fopen(cur_file, 'r');
data = fread(fileID,[1,3600*ADC_SAMPLING_FREQ],'uint16','n');
fclose(fileID);
data(2,:) = bitget(data(1,:),ones(1,length(data))*1,'uint16'); % Newer version injects bit for phase
data(1,:) = bitset(data(1,:),ones(1,length(data))*1,ones(1,length(data))*0,'uint16');
data = transpose(data);

        
if have_next
    disp("loading next overlap");
    fileID = fopen(next_file, 'r');
    data_local = fread(fileID,[1,3600*ADC_SAMPLING_FREQ],'uint16','n');   
    fclose(fileID);
    data_local(2,:) = bitget(data_local(1,:),ones(1,length(data_local))*1,'uint16'); % Newer version injects bit for phase
    data_local(1,:) = bitset(data_local(1,:),ones(1,length(data_local))*1,ones(1,length(data_local))*0,'uint16');
    data = [data; transpose(data_local)];
end

OUTPUT_SAMPLE_RATE = 100;

% Uncalibrated E_field from CU EFM:
E_field = process_hilbert(data, ADC_SAMPLING_FREQ, OUTPUT_SAMPLE_RATE, phase_offset);
timeEFM = datetime(dvec(1),dvec(2),dvec(3),dvec(4),0,0) + seconds( ((0:length(E_field) - 1) )/OUTPUT_SAMPLE_RATE);


%% Apply calibration curve:

E_field_calib = interp1(cal_data.efmVolts, cal_data.E_field_calib, E_field,'linear','extrap');

% figure();
% plot(timeEFM, E_field);
%% Downsample E_field to Campbell rate, for cross-corellation
clf;
fig = figure(2);
hold on;

E_field_down = resample(E_field_calib, 4, OUTPUT_SAMPLE_RATE);
timeEFM_down = timeEFM(1) + seconds( ((0:length(E_field_down) - 1) )/4);
E_field_down(isnan(E_field_down)) = 0;

ts_ind = find(timec(1) == timeEFM_down);
shift = ts_ind - lag;

disp(shift)
disp(ts_ind)


% Fitsky!

leftind = shift;
rightind = shift+length(timec) - 1;

EFD_short =  E_field_down(leftind:rightind);
A = [EFD_short, ones(length(EFD_short),1)];
size(A)
B = EfieldC;
size(B)
fit = A\B;
gain = fit(1);
offset = fit(2);




hold on;
plot(timec, EfieldC);
plot(timec, gain*EFD_short + offset);
legend('Campbell',EFM);
title(sprintf('Cross-calibration Fit at %s',site_name));
ylabel({'Electric field [V/m]';sprintf('Campbell C_{site}=%g',Csite)});
grid on;
fprintf('Gain = %g, offset = %2g\n',gain,offset);

xcoords = get(gca,'xlim');
ycoords = get(gca,'ylim');
xt = xcoords(1) + 0.7*(xcoords(2) - xcoords(1));
yt = ycoords(1) + 0.1*(ycoords(2) - ycoords(1));

text(xt,yt,sprintf('Gain = %2.3f\nOffset = %2.2f',gain,offset));
saveas(gca, fullfile(fig_dir,sprintf('Calibration_curve_%s_%s.png',site_name, EFM)));

