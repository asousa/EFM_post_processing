
a = 120000; b = 121000;
s2 = sig(a:b);
t2 = time(a:b);

figure;
hold on;
for(up_factor=[1,2,5,10,100,200, 500, 1000])
    [su, tu] = resample(s2, t2, up_factor*ADC_SAMPLING_FREQ);
    figure; plot(t2,s2); hold on; plot(tu,su);

    [~, pklocpos] = findpeaks(su);
    [~, pklocneg] = findpeaks(-su);
    [envu,envl] = envPeak(su,pklocpos,pklocneg);
    mag_up = (envu-envl)/2;
%     figure; plot(tu, mag_up);

    [mag, t] = resample(mag_up, tu, ADC_SAMPLING_FREQ);

    [~, pklocpos] = findpeaks(s2);
    [~, pklocneg] = findpeaks(-s2);
    [envu,envl] = envPeak(s2,pklocpos,pklocneg);
    mag_orig = (envu-envl)/2;
%     hold on; plot(t2, mag_orig);
%     plot(t, mag);
end

%% Clean up the phase data if it's noisy!
phase = data((data(:,1)~=0),2);
size(phase)
% p2 = phase(1:end-2).*phase(2:end-1); %.*phase(3:end);
% p2 = (phase(1:end-3) + phase(2:end-2) + phase(3:end-1) + phase(4:end))/4;
phase = square(2*pi*time*ADC_SAMPLING_FREQ); % Generate a synthetic wave
% length(p2)
% p2 = lowpass(phase,100,ADC_SAMPLING_FREQ);
% p2 = round(p2);
% figure; plot(phase(1:100)); hold on; plot(p2(1:100))
% ylim([-.2, 1.2])
% phase = padarray(p2, (length(phase) - length(p2))/2);
length(phase)
% phase = 2.*phase - 1;  % normalize to plus



figure; plot(phase(1:200),'--'); hold on;
plot(p2(1:200));
ylim([-0.2, 1.2]);

%%
nfft = 4096;

spectrogram(sig, blackman(nfft),nfft/2, nfft,ADC_SAMPLING_FREQ);

z = hilbert(sig);
instfrq = ADC_SAMPLING_FREQ/(2*pi)*diff(unwrap(angle(z)));
figure;
plot(time(2:end), instfrq)

%%
[su, tu] = resample(s2, t2, up_factor*ADC_SAMPLING_FREQ);
figure; plot(t2(1:100),s2(1:100)); hold on; plot(tu(1:100*up_factor),su(1:100*up_factor));

%%
map1 = load('CR1000_EFM_arbitrary_voltages Map.mat');
map2 = load('CR1000_EFM_redo Map.mat');

figure; plot(map1.campbell_volts, map1.campbell_field, map2.campbell_table_volts, map2.campbell_field)
