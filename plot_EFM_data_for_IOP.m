% Plot data from a set of EFMs, for an IOP:
close all; clear all;

processed_data_dir= "/Volumes/lairdata/EFM/RELAMPAGO Data/Austin Reprocessed Data/Uncalibrated";
fig_dir = fullfile(processed_data_dir,"IOP figures");
site_names = ["Cordoba","Manfredi","Pilar","Villa-del-Rosario","Villa-Carlos-Paz"];

OUTPUT_SAMPLE_RATE = 100; 



%% time spans of interest (The campaign IOPs)
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
%% Load data for a given span
% IOP_ind = 1;
% for IOP_ind=1:length(spans)
for IOP_ind=9:9
    span = spans(IOP_ind,:);
     % Load all data for a timespan and plot it
    start_date = span(1); end_date = span(2);

    % one file per hour
    dates_to_do = start_date + hours(0:hours(end_date - start_date));
    
    fprintf("Doing span #%d (%s - %s)\n",IOP_ind, start_date, end_date);

    Edata = containers.Map;
    Tdata = containers.Map;

    for site_ind=1:length(site_names)
        site_name = site_names(site_ind);

        Evec = nan(length(dates_to_do)*OUTPUT_SAMPLE_RATE*60*60,1);
        tvec = (0:length(dates_to_do)*60*60*OUTPUT_SAMPLE_RATE - 1)/(60*60*OUTPUT_SAMPLE_RATE);
        for i=1:length(dates_to_do)
            dvec = datevec(dates_to_do(i));
            name = sprintf("%02d.mat",dvec(4));
            disp(name);
            odir = fullfile(processed_data_dir,site_name,sprintf('%d',dvec(1)),sprintf('%d',dvec(2)), sprintf('%d',dvec(3)));

            if isfile(fullfile(odir,name))
                data = load(fullfile(odir,name));
                hr = hours(dates_to_do(i) - dates_to_do(1));

                t_start = hr*60*60*OUTPUT_SAMPLE_RATE + 1;
                t_end = t_start + 60*60*OUTPUT_SAMPLE_RATE - 1;

                Evec(t_start:t_end) = data.E_field_raw;
            end
        end

        Edata(site_name) = Evec;
        Tdata(site_name) = tvec;

    end
    %% Plot it
    close all;
    figure(1);
    axs = [];
    for site_ind=1:length(site_names)
        t_data_labeled = start_date + hours(Tdata(site_name));
        site_name = site_names(site_ind);
        ax = subplot(length(site_names), 1, site_ind);
        plot(ax,t_data_labeled, Edata(site_name),'LineWidth',1);
        ylabel(site_name);
        axs = [axs; ax];
        ylim([-1,1]); % Plotting unscaled instrument units
    end
    linkaxes(axs,'x');

    for i=1:length(axs)-1
    %     axs(i);
        set(axs(i),'XTickLabels',[]);
    end
    xlim([start_date, end_date]);
    sgtitle(sprintf("Uncalibrated EFM data for IOP\n%s - %s",start_date, end_date));
    figname = fullfile(fig_dir,sprintf('IOP_%d_%s-%s.png',IOP_ind,...
        datestr(start_date,'ddmmyy_HH'),...
        datestr(end_date,'ddmmyy_HH')));

%     saveas(gca, figname);
end
