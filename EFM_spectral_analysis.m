% Plot spectra of processed data for each of the times of interest:

root_dir = '/Volumes/lairdata/EFM/RELAMPAGO Data/Campaign Data Processed';
sitenames = ["Cordoba", "Manfredi", "Pilar","VCP", "VDR"];

%%
% time spans of interest
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
% Load data between span for a particular site:
span = spans(1,:)

Edata = containers.Map;
Tdata = containers.Map;

for site_ind=1:length(sitenames)
% for site_ind=1
    site=sitenames(site_ind)
    file_dir  = fullfile(root_dir, site,'*.mat');
    filelist  = dir(file_dir);
    filelist  = filelist(~[filelist.isdir]);  %remove folders from list
    filetimes = datetime.empty(0,1);
    
    for ft_ind=1:length(filelist)
       name = filelist(ft_ind).name;
       try
           dt = datetime(name(1:length(name)-4),'InputFormat','yy-MM-dd','PivotYear',2000);
           filetimes = [filetimes; dt];
       catch
           continue
       end
    end
    filetimes = sort(filetimes);    
    files_to_load = filetimes( (filetimes>=(span(1) - days(1))) & filetimes <= span(2))

    Edata(site) = double.empty(0);
    Tdata(site) = datetime.empty(0);
    
    for i=1:length(files_to_load)
        fpath = fullfile(root_dir, site,strcat(datestr(datetime(files_to_load(i)),'yy-mm-dd', 2000), ".mat"));
        disp('loading ' + fpath);
        d = load(fpath);

        Edata(site) = [Edata(site); (d.E_field_calib)];
        Tdata(site) = [Tdata(site); (d.time)];

    end
end

%% plot it!
ax_list = []
fig = figure(1)
for i=1:2 %length(sitenames)
    ax = subplot(2,1,i)
    site = sitenames(i)
    plot(ax, Tdata(site), Edata(site));
    ax_list = [ax_list; ax];
end
