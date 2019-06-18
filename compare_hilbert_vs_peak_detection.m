close all; clear all;

ADC_SAMPLING_FREQ = 1000;
ADC_REF = 1.81;
up_factor = 10;
% Fiddle with this if we're having issues with large spikes from equal positive to negative values
phase_offset = 1; 
raw_data_dir = "/Volumes/lairdata/EFM/RELAMPAGO Data/Campaign Data";
out_data_dir = "/Volumes/lairdata/EFM/RELAMPAGO Data/Austin Reprocessed Data/Without Site Corrections";
site_name = "Cordoba";

cal_dir = "/Volumes/lairdata/EFM/Field Mill Post-Campaign Calibration/EFM calibration maps 5-12-2019";
% Redo this for each site pls
cal_file = fullfile(cal_dir,"EFM008_map_2019-05-13.mat");  
cal_data = load(cal_file); % returns efmVolts, E_field_calib


%%

%%
input_filepath = "/Volumes/lairdata/EFM/Field Mill Post-Campaign Calibration/Calibration Data 5-12-2019";

fname_full = fullfile(input_filepath, 'EFM008.bin');

disp("Loading raw file");
fileID = fopen(fname_full, 'r');
data = fread(fileID,[1,3600*ADC_SAMPLING_FREQ],'uint16','n'); %% Change 2 to 1 in newer version
fclose(fileID);
data(2,:) = bitget(data(1,:),ones(1,length(data))*1,'uint16'); % Newer version injects bit for phase
data(1,:) = bitset(data(1,:),ones(1,length(data))*1,ones(1,length(data))*0,'uint16');
data = transpose(data);

time = 1/(ADC_SAMPLING_FREQ)*(0:3600*(ADC_SAMPLING_FREQ)-1);
time = transpose(time);


disp("detrending");
%     time = time(data(:,1)~=0);
%     sig = data((data(:,1)~=0),1)/65535*ADC_REF;
sig = data(:,1)/65535*ADC_REF;
% sig = sig - mean(sig);
% sig = detrend(sig,'linear',1:1e7:length(sig));

% sig = lowpass(sig,110,ADC_SAMPLING_FREQ);
sig = bandpass(sig, [99,101], ADC_SAMPLING_FREQ);
sig(data(:,1)==0) = NaN;

disp("upsampling");
[su, tu] = resample(sig, time, up_factor*ADC_SAMPLING_FREQ);

[~, pklocpos] = findpeaks(su);
[~, pklocneg] = findpeaks(-su);
[envu,envl] = envPeak(su,pklocpos,pklocneg);
mag_up = (envu-envl)/2;
[mag, t_down] = resample(mag_up, tu, ADC_SAMPLING_FREQ);

%     E_field = mag(data(:,1)~=0); % No polarity yet
E_field = mag; % No polarity yet
E_field(data(:,1)==0)= NaN;
% Extract Phase Data
%     phase = data((data(:,1)~=0),2);
phase = data(:,2);
phase(data(:,1)==0) = NaN;
phase = 2.*phase - 1;  % normalize to plus/minus 1

% Roll the phase left or right, if needed
phase = circshift(phase,phase_offset); 

% upsample the phase vector as well
[pu, tu] = resample(phase, time, up_factor*ADC_SAMPLING_FREQ);

t_pkloc = [tu(pklocpos);tu(pklocneg)];
phase_pkloc = [(pu(pklocpos)>=0)*1+(pu(pklocpos)<0)*(-1); ...
    (-pu(pklocneg)>=0)*1+(-pu(pklocneg)<0)*(-1)];
pol = interp1(t_pkloc,phase_pkloc,time,'linear','extrap');
pol = sign(pol);

% Apply Polarity and Calibration
E_field = E_field.*(pol);

disp(length(E_field))
disp(length(time))
disp(sum(isnan(E_field)));
% Apply calibration
E_field_calib = interp1(cal_data.efmVolts, cal_data.E_field_calib, E_field,'linear','extrap');
% [folder, baseFileName, extension] = fileparts(flist(i).name);
% outfile = fullfile(odir, baseFileName + ".mat");
% save(outfile,'E_field_calib');

%% 
% % Try it with the hilbert method!
sig = data(:,1)/65535*ADC_REF;
sig = bandpass(sig, [99,101], ADC_SAMPLING_FREQ);
mag = abs(hilbert(sig));
% spectrogram(mag, hanning(nfft),nfft/2,[],ADC_SAMPLING_FREQ,'yaxis');

mag(data(:,1)==0) = NaN;

phase_offset = 2;
phase = data(:,2);
phase = 2.*phase - 1;  % normalize to plus/minus 1

% Roll the phase left or right, if needed
phase = circshift(phase,phase_offset);

phase = bandpass(phase, [99,101], ADC_SAMPLING_FREQ);
% Compute the relative phase between the signal and the encoder
pol_hil = angle(hilbert(sig)) - angle(hilbert(phase));
pol_hil = abs(mod(pol_hil,2*pi) - pi);  % Wrap any angles, scale down by pi
pol_hil = 2*(pol_hil>pi/2) - 1;         % Threshold it at pi/2

figure(2);
% plot(abs(hp))
plot(time, pol_hil)

figure(1);
ax1 = subplot(3,1,1);
plot(ax1, time, E_field);
ax2 = subplot(3,1,2);
plot(ax2, time, mag);
ax3=  subplot(3,1,3);
plot(ax3, time, E_field,'r',time, pol_hil.*mag,'b');

%%
figure(3);
ax1 = subplot(2,1,1);
plot(ax1, time, angle(hilbert(sig))); hold on; plot(time, angle(hilbert(phase)));
ylabel('Negative signal');
legend('Signal phase','Encoder phase');
xlim([150,150.1])
ax2 = subplot(2,1,2);
plot(ax2,time, angle(hilbert(sig))); hold on; plot(time, angle(hilbert(phase)));
ylabel('positive signal');
xlim([1020,1020.1])
%%
figure(4);
ax1 = subplot(2,1,1);
plot(ax1, time, E_field);
ax2 = subplot(2,1,2);
p1 = angle(hilbert(sig)); p2 = angle(hilbert(phase));
% plot(ax2, time, p1,'r',time, p2,'b');
pol_hil = abs(mod(angle(hilbert(sig)) - angle(hilbert(phase)),2*pi) - pi);
pol_hil = 2*(pol_hil>pi/2) - 1;
plot(ax2, time,pol_hil);
linkaxes([ax1, ax2],'x');
% ylim([-360,360]);
% xlim([150, 150.1]);



