addpath(genpath('/Users/barrab01/Documents/Repos/matnwb'));
addpath(genpath('/Users/barrab01/Documents/Repos/ephysMH'));
addpath(genpath('/Users/barrab01/Documents/Repos/chronux_2_12'));
clc
clear all

folder = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_14/';
filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_14/040108__24_11_14.nwb';
% 
folder ='/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_16/';
filename ='/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_16/040108__24_11_16.nwb';
%  
folder = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040116/24_11_13/';
filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040116/24_11_13/040116__24_11_13.nwb';
% 
folder = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/12_11_15/';
filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/12_11_15/040118__24_11_15.nwb';
% 
folder = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/24_11_18/';
filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/24_11_18/040118__24_11_18.nwb';
% 
folder = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_11/';
filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_11/040136__24_11_11.nwb';
% 
folder = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_12/'; 
filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_12/040136__24_11_12.nwb'; 

savefilename = fullfile(folder, '/reduced_dataset.mat');
savedistrialname = fullfile(folder, 'discarded_trials.mat');
savetrialnumbers = fullfile(folder, 'trials_numbers.mat');
savefigurename = fullfile(folder, 'PC_qualitycheck.png');
nwbFile = nwbRead(filename);
analogFS = 1000; 
dt = 1/analogFS; 
odorlabels = {'ethyltiglate', 'ethylbutyrate', 'acetophenone', 'heptanone'}; 


%% Read data 

% Retrieve start and stop times
trials_start_stop_time = [nwbFile.intervals.get('odor_trials').start_time.data(:), ...
                        nwbFile.intervals.get('odor_trials').stop_time.data(:)];

% Extract trial
Odors = nwbFile.intervals.get('odor_trials').vectordata.get('pins').data(:); 
Concentration = nwbFile.intervals.get('odor_trials').vectordata.get('concentration').data(:); 
FVO = nwbFile.intervals.get('odor_trials').vectordata.get('fvo').data(:); 
FVC = nwbFile.intervals.get('odor_trials').vectordata.get('fvc').data(:); 
% Analog Sniff signals
Sniff = nwbFile.acquisition.get('respiration_data').data(:);  
SniffTime = nwbFile.acquisition.get('respiration_data').timestamps(:);  
Sniffdt = mean(diff(nwbFile.acquisition.get('respiration_data').timestamps(:))); 
fSniff = lowPass(Sniff, round(1/Sniffdt), 10, 3); 
baselineSniff = lowPass(Sniff, round(1/Sniffdt),0.8, 2); 
SniffFS = 1/mean(diff(nwbFile.acquisition.get('respiration_data').timestamps(:))); 
% Sniff features
Prex = nwbFile.analysis.get('respiratory_features').vectordata.get('prex').data(:); 
Postx = nwbFile.analysis.get('respiratory_features').vectordata.get('postx').data(:); 

figure
times = [1, 2]
idx_plot = floor(times.*SniffFS); 
[temp, inh_in_window] = arrayfun(@(t) min(abs(Prex-t)), times); 
plot(SniffTime(idx_plot(1): idx_plot(2)), Sniff(idx_plot(1): idx_plot(2)))
hold on

arrayfun(@(t) xline(t, '--r', 'LineWidth', 1.5), Prex(inh_in_window(1):inh_in_window(2)));
arrayfun(@(t) xline(t, '--r', 'LineWidth', 1.5), Postx(inh_in_window(1):inh_in_window(2)));


%% Units quality check
%Define each QC criterion separately
location_drop = 'not in';  % Example: location filter for specific areas
ks_label = 'good';  % Example: Kilosort quality
SI_label = 'good';  % Example: SI quality
cutoff_hz = 0.2;  % Example: firing rate cutoff
cutoff_amp = 10 ; % Example: amplitude cutoff
max_above = 400 ; % Example: max depth for 'PC:' locations
% Compute depth
temp = cellfun(@(x) split(x, ':'), nwbFile.units.vectordata.get('location').data(:), 'UniformOutput', false); 
depth = cellfun(@(x) str2num(x{2}), temp); 

% General quality check (same for PC and OB)
quality_units = (cellfun(@(x) ~contains( x , location_drop ), nwbFile.units.vectordata.get('location').data(:)) & ...
    cellfun(@(x) contains( x , ks_label ), nwbFile.units.vectordata.get('quality').data(:)) & ...
    cellfun(@(x) contains( x , SI_label ), nwbFile.units.vectordata.get('SI_quality').data(:)) &...
    nwbFile.units.vectordata.get('firing_rate').data(:)>cutoff_hz & ...
    (nwbFile.units.vectordata.get('amplitude').data(:)>cutoff_amp)'); 

PC_units = (cellfun(@(x) contains( x , 'PC' ), nwbFile.units.vectordata.get('location').data(:)) & ...
    depth>max_above & ...
    quality_units);
PC_units_idx = find(PC_units); 
OB_units = (cellfun(@(x) contains( x , 'OB' ), nwbFile.units.vectordata.get('location').data(:)) & ...
    quality_units); 
OB_units_idx = find(OB_units); 
%disp([sum(PC_units), sum(OB_units)])

%% Extract trial numbers per condition
iodor = 0; 
for odor = unique(Odors)'
    iodor = iodor +1; 
    trials_numbers.(odorlabels{iodor}) = struct(); 
    idx_odor = find(Odors ==odor);
    all_conc = unique(Concentration(idx_odor)); 
    trials_numbers.(odorlabels{iodor}).concs = all_conc; 
    trials_numbers.(odorlabels{iodor}).trial_idxs = []; 
    for iconc  = 1: length(all_conc)
        trials_numbers.(odorlabels{iodor}).trial_idxs{iconc} = intersect(find(Concentration == all_conc(iconc)),idx_odor );
         
    end
    
end
%save(savetrialnumbers,"trials_numbers","-v7.3")

%% Visualize respiration traces alignment
o = 8; % acetophenone
c = 0.1;
figure
idx_trials = find(Odors == o & Concentration == c); 
ii = 0;
iii = 0; 
vertical_offset = 0.3; 
for itrial =  idx_trials'
    iii = iii +1; 
    ii = ii -vertical_offset;
    start_time = FVO(itrial) - 10;
    end_time = FVO(itrial) + 6;
    start_idx = searchsorted(SniffTime, start_time);
    end_idx = searchsorted(SniffTime, end_time);
    % Extract respiration traces
    respiration_trace_window = fSniff(start_idx:end_idx)+ii; 
    postfvoinh = find(Prex>FVO(itrial)); 
    inh1 = Prex(postfvoinh(1)); 
    inh2 = Prex(postfvoinh(2)); 

    % Normalize time axis relative to FVO
    trace_times = SniffTime(start_idx:end_idx) - FVO(itrial); 
    
    hold on
    fill([0, FVC(itrial)-FVO(itrial), FVC(itrial)-FVO(itrial), 0 ],...
        [mean(respiration_trace_window)-vertical_offset/2, mean(respiration_trace_window)-vertical_offset/2, ...
        mean(respiration_trace_window)+vertical_offset/2, mean(respiration_trace_window)+vertical_offset/2],  'b', 'FaceAlpha', 0.4, 'EdgeColor', 'b', 'EdgeAlpha', 0.2)

    fill([inh1-FVO(itrial), inh2-FVO(itrial), inh2-FVO(itrial), inh1-FVO(itrial) ],...
        [mean(respiration_trace_window)-vertical_offset/2, mean(respiration_trace_window)-vertical_offset/2, ...
        mean(respiration_trace_window)+vertical_offset/2, mean(respiration_trace_window)+vertical_offset/2],  [0.8, 0.8, 0.8], 'FaceAlpha', 1, 'EdgeColor', [0.8, 0.8, 0.8], 'EdgeAlpha', 0.5)
    
    plot(trace_times, respiration_trace_window,'k', 'Linewidth', 1.5)
    title(['Vial = ', num2str(o), ',  Conc = ', num2str(c)])
    text(trace_times(end) + 0.2, mean(respiration_trace_window), ['Trial # ', num2str(itrial)])
    percent_complete = (iii / length(idx_trials)) * 100;
    disp([' ----   ', num2str(percent_complete), ' % of Sniff plotting run   ----'])
end

%% Count spikes in 300 ms window after inhalation


all_odors = unique(Odors); 
 
spikes = nwbFile.units.spike_times.data(:); 
mywindow = 0.3; 
FR_table = cell(length(PC_units_idx), 4*8 );
baseline_FR_table = cell(length(PC_units_idx), 4*8 );
for i = 1: length(PC_units_idx)
    for j = 1 : 8
        FR_table{i,j} = []; 
        baseline_FR_table{i,j} = []; 
    end 
end
figure
for itrial  = 1 : length(Odors)
    
    % Find odor and concentration indexes to pull together trials with the
    % same conditions
    iodor = find(all_odors == Odors(itrial)); 
    idx_odor = find(Odors ==iodor+4); 
    all_conc_thisodor = unique(Concentration(idx_odor));
    iconc = find(all_conc_thisodor == Concentration(itrial)); 
    condition_idx = (iodor-1)*length(all_conc_thisodor) + iconc; 
    %disp([iodor, iconc, condition_idx]);
    

    
    all_prex_before_fvo_idx = find(Prex<FVO(itrial));
    all_prex_after_fvo_idx = find(Prex>FVO(itrial));
    all_postx_before_fvo_idx=  find(Postx<FVO(itrial));
    % Start and end of baseline
    start_bas = Prex(all_prex_before_fvo_idx(end-2));  
    end_bas = start_bas + mywindow; Postx(all_postx_before_fvo_idx(end)); 
    
    % Start and end of inh
    
    start_inh = Prex(all_prex_after_fvo_idx(1)); 
    end_inh =start_inh + mywindow; 
    
    % Plot to check 
%     start_time = FVO(itrial) -1.5;
%     end_time = FVO(itrial) + 1.5;
%     start_idx = searchsorted(SniffTime, start_time);
%     end_idx = searchsorted(SniffTime, end_time)
%     % Extract respiration traces
%     respiration_trace_window = fSniff(start_idx:end_idx); 
%     trace_times = SniffTime(start_idx:end_idx) - FVO(itrial); 
%     plot(trace_times, respiration_trace_window,'k', 'Linewidth', 1.5)
%     hold on
%     plot([start_inh-FVO(itrial), start_inh-FVO(itrial)], [min(respiration_trace_window), max(respiration_trace_window)], 'r--')
%     plot([end_inh-FVO(itrial), end_inh-FVO(itrial)], [min(respiration_trace_window), max(respiration_trace_window)], 'r--')
%     plot([start_bas-FVO(itrial), start_bas-FVO(itrial)], [min(respiration_trace_window), max(respiration_trace_window)], 'b--')
%     plot([end_bas-FVO(itrial), end_bas-FVO(itrial)], [min(respiration_trace_window), max(respiration_trace_window)], 'g--')
%     plot([0, 0], [min(respiration_trace_window), max(respiration_trace_window)], 'k--')
%     [FVO(itrial), start_bas, end_bas, start_inh, end_inh]
%     pause
%     clf
    % Loop through units and find spikes in baseline and first sniff
    iiunit = 0; 
    for iunit = PC_units_idx'
        iiunit = iiunit+1; 
        last_spike = nwbFile.units.spike_times_index.data(iunit); 
        if iunit >1
            first_spike = nwbFile.units.spike_times_index.data(iunit-1)+1; 
        else
            first_spike = 1; 
        end
        allspikes = spikes(first_spike:last_spike); 
        % First sniff 
        FR = length(find( allspikes> start_inh & allspikes< end_inh))/ mywindow; 
        FR_table{iiunit, condition_idx} = [FR_table{iiunit, condition_idx}, FR]; 
        % Baseline
        baseline_FR = length(find( allspikes> start_bas & allspikes< end_bas))/ mywindow; 
        baseline_FR_table{iiunit, condition_idx} = [baseline_FR_table{iiunit, condition_idx}, baseline_FR]; 
    end
    percent_complete = (itrial / length(Odors)) * 100;
    disp([' ----   ', num2str(percent_complete), ' % of spike extraction and heatmaps computations run   ----'])

end
% Average FR table 
FRMAT = cellfun(@(x) mean(x), FR_table); 
bas_FRMAT = cellfun(@(x) mean(x), baseline_FR_table);
mean_baseline_byunit = mean(bas_FRMAT, 2);
std_baseline_byunit = std(bas_FRMAT')'; 
mean_activity_byunit = mean(FRMAT, 2);
valid_units = find(std_baseline_byunit~=0); 

% Loop across stimulus condition again  and units and zscore each unit (not
% do it for all trials actually
disp('Computing z-score ...')
z_scored_spike_rates = {}; 
for icond = 1 : size(FR_table, 2)
    for iunit = valid_units'
        for irep = 1 : length(FR_table(iunit, icond))
            z_scored_spike_rates{iunit, icond}(irep) = (FR_table{iunit, icond}(irep) - mean_baseline_byunit(iunit)) / std_baseline_byunit(iunit);
        end
    end
end
ZFRMAT = cellfun(@(x) mean(x), z_scored_spike_rates); 
disp('Done')
disp('Computing response index...')
% Compute response index
figure
for icond = 1 : size(FR_table, 2)
    for iunit = valid_units'
        
        baseline_rates = baseline_FR_table{iunit, icond}; 
        stimulus_rates = FR_table{iunit, icond}; 
        labels = [repmat([0], 1, length(baseline_rates)) ,  repmat([1], 1,length(stimulus_rates))];  
        scores = [baseline_rates , stimulus_rates];    
        %disp([labels; scores])
        [X,Y,T,AUC] = perfcurve(labels,scores, 1); 
        plot(X, Y)
        hold on
        response_index = 2 * (AUC - 0.5);
        response_index_MAT(iunit, icond) = response_index; 
    end
    %percent_complete = (itrial / length(Odors)) * 100;
    %disp([' ----   ', num2str(percent_complete), ' % of response index computations run   ----'])
    pause(0.01)
end
disp('Done')

%% Plotting heatmaps
disp('Saving figure...')
blueredcmap = blue_white_red(100); %flipud(colormap('summer')); 
bluecmap =flipud(blueredcmap(1:50,:));  %flipud(colormap('bone')); 

fig = figure; 
%imagesc(bas_FRMAT)
ax(1) = subplot(1, 3, 1);
imagesc(FRMAT)
colormap(ax(1),bluecmap);
colorbar();
hold on
xline([8.5,16.5,24.5,32.5], 'k--', 'Linewidth', 2)
title('Firing rates')
xlabel('Stim Conditions')
ylabel('#Unit')
ax(2) = subplot(1, 3, 2);
minv = quantile(ZFRMAT(:), [0.05]); 
maxv = quantile(ZFRMAT(:), [0.95]); 
maxV = max([minv, maxv]); 
imagesc(ZFRMAT, [-maxV, maxv]);
colormap(ax(2),blueredcmap);
colorbar();
hold on
xline([8.5,16.5,24.5,32.5], 'k--', 'Linewidth', 2);
title('Z-scored firing rate');
xlabel('Stim Conditions');
ylabel('#Unit');

ax(3) = subplot(1, 3, 3);
%minv = quantile(response_index_MAT(:), [0.05]); 
%maxv = quantile(response_index_MAT(:), [0.95]); 
imagesc(response_index_MAT); %[minv, maxv]
colormap(ax(3),blueredcmap);
colorbar();
hold on
xline([8.5,16.5,24.5,32.5], 'k--', 'Linewidth', 2);
title('Response Index');
xlabel('Stim Conditions');
ylabel('#Unit');
saveas(fig, savefigurename);
disp('Done')

