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
cal_file = fullfile(cal_dir,"EFM011_map_2019-05-13.mat");  
cal_data = load(cal_file); % returns efmVolts, E_field_calib


%%

file_date = datetime(2018,12,4,0,0,0 );
dvec = datevec(file_date);
disp(dvec)
fpath = fullfile(raw_data_dir,site_name,'DATA',...
sprintf('%d',dvec(1)),sprintf('%d',dvec(2)), sprintf('%d',dvec(3)),'*.bin');
disp(fpath)
flist = dir(fpath);


%% Prepare output directory
odir = fullfile(out_data_dir,site_name,sprintf('%d',dvec(1)),sprintf('%d',dvec(2)), sprintf('%d',dvec(3)));
mkdir(odir);


%%

for i=1:1 %length(flist)
    disp(flist(i).name);
    fname_full = fullfile(raw_data_dir,site_name,'DATA',...
                sprintf('%d',dvec(1)),sprintf('%d',dvec(2)), sprintf('%d',dvec(3)),...
                flist(i).name);
    
    fileID = fopen(fname_full, 'r');
    data = fread(fileID,[1,3600*ADC_SAMPLING_FREQ],'uint16','n'); %% Change 2 to 1 in newer version
    fclose(fileID);
    data(2,:) = bitget(data(1,:),ones(1,length(data))*1,'uint16'); % Newer version injects bit for phase
    data(1,:) = bitset(data(1,:),ones(1,length(data))*1,ones(1,length(data))*0,'uint16');
    data = transpose(data);
    
    time = 1/(ADC_SAMPLING_FREQ)*(0:3600*(ADC_SAMPLING_FREQ)-1);
    time = transpose(time);


    
%     time = time(data(:,1)~=0);
%     sig = data((data(:,1)~=0),1)/65535*ADC_REF;
    sig = data(:,1)/65535*ADC_REF;
    sig = detrend(sig,'linear',1:1e7:length(sig));

    sig = lowpass(sig,110,ADC_SAMPLING_FREQ);
    sig(data(:,1)==0) = NaN;

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
    disp(flist(i).name)
    disp(length(E_field))
    disp(length(time))
    disp(sum(isnan(E_field)));
    % Apply calibration
    E_field_calib = interp1(cal_data.efmVolts, cal_data.E_field_calib, E_field,'linear','extrap');
    [folder, baseFileName, extension] = fileparts(flist(i).name);
    outfile = fullfile(odir, baseFileName + ".mat");
%     save(outfile,'E_field_calib');
    
end

%%
flist = dir(fullfile(odir,"*.mat"));
disp(flist)

Evec = zeros(length(flist)*ADC_SAMPLING_FREQ*60*60,1);
tvec = (0:length(flist)*60*60*ADC_SAMPLING_FREQ - 1)/(60*60*ADC_SAMPLING_FREQ);
for i=1:length(flist)
    disp(flist(i).name);
    data = load(fullfile(odir,flist(i).name));
    [folder, baseFileName, extension] = fileparts(flist(i).name);
    hr = str2num(baseFileName);
    t_start = hr*60*60*ADC_SAMPLING_FREQ + 1
    t_end = t_start + 60*60*ADC_SAMPLING_FREQ - 1
    
    Evec(t_start:t_end) = data.E_field_calib;
%     Evec = [Evec; data.E_field_calib];
end

%% This version to load all the files for a single day, and then process it

file_date = datetime(2018,12,4,0,0,0 );
dvec = datevec(file_date);
disp(dvec)
fpath = fullfile(raw_data_dir,site_name,'DATA',...
sprintf('%d',dvec(1)),sprintf('%d',dvec(2)), sprintf('%d',dvec(3)),'*.bin');
disp(fpath)
flist = dir(fpath);


data = zeros(length(flist)*60*60*ADC_SAMPLING_FREQ,2);
time = zeros(length(flist)*60*60*ADC_SAMPLING_FREQ,1);
for i=1:length(flist)
    disp(flist(i).name);
    fname_full = fullfile(raw_data_dir,site_name,'DATA',...
                sprintf('%d',dvec(1)),sprintf('%d',dvec(2)), sprintf('%d',dvec(3)),...
                flist(i).name);
    
    fileID = fopen(fname_full, 'r');
    data_local = fread(fileID,[1,3600*ADC_SAMPLING_FREQ],'uint16','n'); %% Change 2 to 1 in newer version
    fclose(fileID);
    data_local(2,:) = bitget(data_local(1,:),ones(1,length(data_local))*1,'uint16'); % Newer version injects bit for phase
    data_local(1,:) = bitset(data_local(1,:),ones(1,length(data_local))*1,ones(1,length(data_local))*0,'uint16');
    data_local = transpose(data_local);
    
    [folder, baseFileName, extension] = fileparts(flist(i).name);
    hr = str2num(baseFileName);
    
    time_local = hr*60*60 + 1/(ADC_SAMPLING_FREQ)*(0:3600*(ADC_SAMPLING_FREQ)-1);
    time_local = transpose(time_local);

    t_start = hr*60*60*ADC_SAMPLING_FREQ + 1
    t_end = t_start + 60*60*ADC_SAMPLING_FREQ - 1
    
    data(t_start:t_end,:) = data_local;
    time(t_start:t_end) = time_local;
end

%%


sig = data(:,1)/65535*ADC_REF;
% sig = sig - mean(sig);
% sig = detrend(sig,'linear',1:1e7:length(sig));

% sig = lowpass(sig,110,ADC_SAMPLING_FREQ);
sig = bandpass(sig, [99,101], ADC_SAMPLING_FREQ);
sig(data(:,1)==0) = NaN;

[su, tu] = resample(sig, time, up_factor*ADC_SAMPLING_FREQ);
figure(1);
t1 = 1000;
t2 = 5000;

plot(tu(t1:t2),su(t1:t2));

%%
nfft = 1024;

t1 = 1;
t2 = length(sig);
[s,f,t] = spectrogram(sig(t1:t2), blackman(nfft),nfft/2, [], ADC_SAMPLING_FREQ);

% find peak frequencies:
[a, pf_loc] = max(log10(abs(s)));
% disp(f(pf_loc))

% sig_demod = zeros(length(sig),1);
% 
% for i=1:length(t)
%     t1 = (nfft/2)*(i-1) + 1;
%     t2 = (nfft/2)*i;
%     invec = sig(t1:t2);
%     fmax = f(pf_loc(i));
%     disp(fmax);
%     sig_demod(t1:t2) = amdemod(invec, fmax, ADC_SAMPLING_FREQ);
% end

figure(1);
pcolor(t,f,10*log10(abs(s)));
shading flat;

sh = hilbert(su);
[sig_demod, t_down] = resample(abs(sh), tu, ADC_SAMPLING_FREQ);

sig2 = abs(hilbert(sig));

figure(2);
% % plot(time, sig_demod);
% ax1 = subplot(2,1,1);
% plot(ax1,time, mag);
% ax2 = subplot(2,1,2);
% plot(ax2,time, sig_demod);

plot(time, mag, 'r',time,sig_demod,'b', time, sig2,'c');
