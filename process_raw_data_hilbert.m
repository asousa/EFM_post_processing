close all; clear all;

ADC_SAMPLING_FREQ = 1000;
OUTPUT_SAMPLE_RATE = 100; % Hz, what we're decimating to

% Fiddle with this if we're having issues with large spikes from equal positive to negative values
phase_offset = 1; 
raw_data_dir = "/Volumes/lairdata/EFM/RELAMPAGO Data/Campaign Data";
out_data_dir = "/Volumes/lairdata/EFM/RELAMPAGO Data/Austin Reprocessed Data/Uncalibrated";
% cal_dir = "/Volumes/lairdata/EFM/Field Mill Post-Campaign Calibration/EFM calibration maps 5-12-2019";

% site_name = "Cordoba"; % cal_file = fullfile(cal_dir,"EFM011_map_2019-05-13.mat");  
% site_name = "Manfredi"; %cal_file = fullfile(cal_dir,"EFM004_map_2019-05-13.mat");  
% site_name = "Pilar";   % cal_file = fullfile(cal_dir,"EFM006_map_2019-05-13.mat");  
% site_name = "Villa-del-Rosario"; %cal_file = fullfile(cal_dir,"EFM002_map_2019-05-13.mat");  
site_name = "Villa-Carlos-Paz";  %cal_file = fullfile(cal_dir,"EFM008_map_2019-05-13.mat");  

% cal_data = load(cal_file); % returns efmVolts, E_field_calib

%%
% time spans of interest (The campaign IOPs)
spans = datetime.empty(0,2);
spans = [spans; [datetime(2018,12,4,11,0,0 ), datetime(2018,12,5,10,0,0)] ];
spans = [spans; [datetime(2018,11,3,13,0,0 ), datetime(2018,11,4,11,0,0)] ];
spans = [spans; [datetime(2018,12,11,16,0,0), datetime(2018,12,11,22,0,0)] ];
spans = [spans; [datetime(2018,11,25,20,0,0), datetime(2018,11,27,20,0,0)] ];
spans = [spans; [datetime(2018,11,4,20,0,0 ), datetime(2018,11,7,10,0,0 )] ];
spans = [spans; [datetime(2018,11,29,14,0,0), datetime(2018,12,1,10,0,0 )] ];
spans = [spans; [datetime(2018,11,21,22,0,0), datetime(2018,11,22,23,0,0)] ];
spans = [spans; [datetime(2018,12,5,15,0,0 ), datetime(2018,12,6,4,0,0  )] ];
spans = [spans; [datetime(2018,12,13,16,0,0), datetime(2018,12,14,8,0,0 )] ];
spans = [spans; [datetime(2018,11,10,15,0,0), datetime(2018,11,13,6,0,0 )] ];
spans = [spans; [datetime(2018,11,2,23,0,0 ), datetime(2018,11,3,2,0,0  )] ];

%%

for s_ind=1:length(spans)
    start_date = spans(s_ind,1);
    end_date = spans(s_ind,2);
    fprintf("Span %d\n",s_ind);


    % start_date = datetime(2018,12,4,0,0,0 );
    % end_date = datetime(2018,12,6,0,0,0 );

    % one file per hour
    dates_to_do = start_date + hours(0:hours(end_date - start_date));
    overlap_samples = 5*ADC_SAMPLING_FREQ; % How many samples to overlap from the adjacent files

    for i=1:length(dates_to_do)
        % Load the main file, with some overlap at the beginning and end from
        % the adjacent files (Hilbert transform takes a bit to settle)
        file_time = dates_to_do(i);
        prev_file_time = file_time - hours(1);
        next_file_time = file_time + hours(1);

        dvec = datevec(prev_file_time);
        prev_file = fullfile(raw_data_dir,site_name,'DATA',...
                    sprintf('%d',dvec(1)),sprintf('%d',dvec(2)), sprintf('%d',dvec(3)),...
                    sprintf('%02d.bin',dvec(4)));
        dvec = datevec(next_file_time);
        next_file = fullfile(raw_data_dir,site_name,'DATA',...
                    sprintf('%d',dvec(1)),sprintf('%d',dvec(2)), sprintf('%d',dvec(3)),...
                    sprintf('%02d.bin',dvec(4)));
        dvec = datevec(file_time);
        cur_file = fullfile(raw_data_dir,site_name,'DATA',...
                    sprintf('%d',dvec(1)),sprintf('%d',dvec(2)), sprintf('%d',dvec(3)),...
                    sprintf('%02d.bin',dvec(4)));
        
        if ~isfile(cur_file)
            fprintf("No data available for %s at %s\n",site_name, file_time);
            continue;
        end
        have_prev   = isfile(prev_file);
        have_next   = isfile(next_file);

        % Load current file
        disp("Loading current file");
        fileID = fopen(cur_file, 'r');
        data = fread(fileID,[1,3600*ADC_SAMPLING_FREQ],'uint16','n');
        fclose(fileID);
        data(2,:) = bitget(data(1,:),ones(1,length(data))*1,'uint16'); % Newer version injects bit for phase
        data(1,:) = bitset(data(1,:),ones(1,length(data))*1,ones(1,length(data))*0,'uint16');
        data = transpose(data);

        % Load previous overlap
        if have_prev
            disp("loading previous overlap");
            fileID = fopen(prev_file, 'r');
            fseek(fileID, overlap_samples, 1);
            data_local = fread(fileID,[1,overlap_samples],'uint16','n');
            fclose(fileID);
            data_local(2,:) = bitget(data_local(1,:),ones(1,length(data_local))*1,'uint16'); % Newer version injects bit for phase
            data_local(1,:) = bitset(data_local(1,:),ones(1,length(data_local))*1,ones(1,length(data_local))*0,'uint16');
            data = [transpose(data_local); data];
        end

        % Load next overlap
        if have_next
            disp("loading next overlap");
            fileID = fopen(next_file, 'r');
            data_local = fread(fileID,[1,overlap_samples],'uint16','n');
            fclose(fileID);
            data_local(2,:) = bitget(data_local(1,:),ones(1,length(data_local))*1,'uint16'); % Newer version injects bit for phase
            data_local(1,:) = bitset(data_local(1,:),ones(1,length(data_local))*1,ones(1,length(data_local))*0,'uint16');
            data = [data; transpose(data_local)];
        end


    %     % Decode the data (moved into a function so we can reuse it elsewhere)
    %     E_field_calib = process_hilbert(data, ADC_SAMPLING_FREQ, ADC_REF,...
    %                                     EFM_local_gain, cal_data, OUTPUT_SAMPLE_RATE, phase_offset);
        E_field_raw = process_hilbert(data, ADC_SAMPLING_FREQ, OUTPUT_SAMPLE_RATE, phase_offset);

        % Trim off any extra from the overlap:
        trim_length = (overlap_samples*OUTPUT_SAMPLE_RATE/ADC_SAMPLING_FREQ);
        if have_prev
            disp("trimming prev");
            E_field_raw = E_field_raw(trim_length + 1:end);
        end
        if have_next
            disp("trimming next");
           E_field_raw = E_field_raw(1:end - trim_length); 
        end

        % Chat about it
        fprintf("File %s: \n",cur_file);
        fprintf("Length(output)=%d\n",length(E_field_raw));
        fprintf("NaNs in input: %d, NaNs in output: %d \n",sum(data(:,1)==0),sum(isnan(E_field_raw)))
        fprintf("Total dropouts: %3.2g seconds \n",(sum(isnan(E_field_raw))/OUTPUT_SAMPLE_RATE))



        odir = fullfile(out_data_dir,site_name,sprintf('%d',dvec(1)),sprintf('%d',dvec(2)), sprintf('%d',dvec(3)));
        if ~isfolder(odir)
            sprintf("Making directory %s\n",odir);
            mkdir(odir);
        end

        outfile = fullfile(odir, sprintf("%02d.mat",dvec(4)));
        save(outfile,'E_field_raw');
    end

end

 %% Load all data for a timespan and plot it

% one file per hour
dates_to_do = start_date + hours(0:hours(end_date - start_date));

dvec = datevec(start_date);
odir = fullfile(out_data_dir,site_name,sprintf('%d',dvec(1)),sprintf('%d',dvec(2)), sprintf('%d',dvec(3)));
% flist = dir(fullfile(odir,"*.mat"));
% disp(flist)

Evec = zeros(length(dates_to_do)*OUTPUT_SAMPLE_RATE*60*60,1);
tvec = (0:length(dates_to_do)*60*60*OUTPUT_SAMPLE_RATE - 1)/(60*60*OUTPUT_SAMPLE_RATE);
for i=1:length(dates_to_do)
    dvec = datevec(dates_to_do(i));
    name = sprintf("%02d.mat",dvec(4));
    disp(name);
    odir = fullfile(out_data_dir,site_name,sprintf('%d',dvec(1)),sprintf('%d',dvec(2)), sprintf('%d',dvec(3)));
    data = load(fullfile(odir,name));
    hr = hours(dates_to_do(i) - dates_to_do(1));
    
%     [folder, baseFileName, extension] = fileparts(flist(i).name);
%     hr = str2double(baseFileName);
    t_start = hr*60*60*OUTPUT_SAMPLE_RATE + 1;
    t_end = t_start + 60*60*OUTPUT_SAMPLE_RATE - 1;
    
    Evec(t_start:t_end) = data.E_field_raw;
%     Evec = [Evec; data.E_field_calib];
end

figure(1);
plot(tvec, Evec)
hold on;
plot(tvec(isnan(Evec)), ones(sum(isnan(Evec)),1)*0,'co');
