%  Bea's paths
addpath(genpath('/Users/barrab01/Documents/Repos/PiriformAnalysis/Dependencies'))
addpath(genpath('/Users/barrab01/Documents/Repos/chronux_2_12'))
addpath(genpath('/Users/barrab01/Documents/Repos/fieldtrip'))
addpath(genpath('/Users/barrab01/Documents/Repos/matnwb'))
rmpath('/Users/barrab01/Documents/Repos/fieldtrip/external/signal/')
%% Load data

mousefolder = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_14/'; 

% Read start and stop of trials
filename_units = [mousefolder, '040108__24_11_14.nwb']; 
nwbFile = nwbRead(filename_units);

% Read important variables
trials_start_stop_time = [nwbFile.intervals.get('odor_trials').start_time.data(:), ...
                        nwbFile.intervals.get('odor_trials').stop_time.data(:)]
Odors = nwbFile.intervals.get('odor_trials').vectordata.get('pins').data(:); 
Concentration = nwbFile.intervals.get('odor_trials').vectordata.get('concentration').data(:); 

% Analog Sniff signals
Sniff = nwbFile.acquisition.get('respiration_data').data(:);  
SniffTime = nwbFile.acquisition.get('respiration_data').timestamps(:);  
Sniffdt = mean(diff(nwbFile.acquisition.get('respiration_data').timestamps(:))); 
%fSniff = lowPass(Sniff, round(1/Sniffdt), 10, 3); 
%baselineSniff = lowPass(Sniff, round(1/Sniffdt),0.8, 2); 
SniffFS = 1/mean(diff(nwbFile.acquisition.get('respiration_data').timestamps(:))); 
% Sniff features
Prex = nwbFile.analysis.get('respiratory_features').vectordata.get('prex').data(:); 
Postx = nwbFile.analysis.get('respiratory_features').vectordata.get('postx').data(:); 

    
FVO = nwbFile.intervals.get('odor_trials').vectordata.get('fvo').data(:); 
% Load LFP


odorlabels = {'ethyltiglate', 'ethylbutyrate', 'acetophenone', 'heptanone'}; 
freq = logspace(log10(0.5),log10(200),102); 
time = [-0.5, 2]

xticks_time = linspace(time(1), time(2), 5);

% Select one every 10 values for frequency
yticks_freq = freq([1:10:length(freq)]);
yticks_pos = [1:10:length(freq)];

%% plot spectrogram
% List how many files in mousefolder
% Get the list of LFP files to process
lfpFiles = dir(fullfile(mousefolder, 'LFP/processed', 'stats_ch*.mat'));
% Write loop iterating on file and keep ch as index

for ch = 1 : 2%length(lfpFiles)
    filename = [mousefolder, 'LFP/processed/stats_ch', num2str(ch), '.mat']
    load(filename)

    figure
    % odor 1
   
    iodor = 0; 
    for odor = unique(Odors)'
        iodor = iodor+1; 
        idx_odor = find(Odors == odor); 
        iconc = 0; 
        for concentration = unique(Concentration(idx_odor))'
            iconc = iconc+1; 
            % Process data for the current concentration
            idx_conc = find(Concentration == concentration);
            idx_trial = intersect(idx_odor, idx_conc);
            lfp = [];
            ii =0; 
            for idx = idx_trial'
                ii = ii+1; 
                lfp(:,:,ii) = stats.Z{idx};
                
                
            end
            subplot(4, 10, (iodor-1)*10 + iconc)
            imagesc(flipud(mean(lfp, 3)), [0, 2])
            
            
            % put xticks 
            xticks_pos = linspace(1, size(lfp, 2), 5); 
            set(gca, 'XTick', xticks_pos, 'XTickLabel', round(xticks_time, 2));
            set(gca, 'YTick', yticks_pos, 'YTickLabel', round(yticks_freq, 2));
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            title([odorlabels{iodor}, ' conc = ', num2str(concentration) ])
            
        
        end
    end
    clear stats
end

%% PLot overall activity post sniff 
% List how many files in mousefolder
% Get the list of LFP files to process
lfpFiles = dir(fullfile(mousefolder, 'LFP/processed', 'stats_ch*.mat'));
% Write loop iterating on file and keep ch as index
yoffset = 1; 
for ch = 1 : 2%length(lfpFiles)
    filename = [mousefolder, 'LFP/processed/stats_ch', num2str(ch), '.mat']
    load(filename)
    % Look for maximums values across all dataset 
    maxactivity = 0; 
    for i = 1 : length(stats.Z)
        temp = stats.Z{i}; 
        maxactivity = max([maxactivity, max(temp(:))]); 
    end
    
    % odor 1
    total_activity = []
    iodor = 0; 
    for odor = unique(Odors)'
        figure
        iodor = iodor+1; 
        idx_odor = find(Odors == odor); 
        iconc = 0; 
        for concentration = unique(Concentration(idx_odor))'
            iconc = iconc+1; 
            % Process data for the current concentration
            idx_conc = find(Concentration == concentration);
            idx_trial = intersect(idx_odor, idx_conc);
            lfp = [];
            ii =0; 
            
            for idx = idx_trial'
                ii = ii+1; 
                total_activity(ii, : ) = mean(stats.Z{idx}, 1); 
                subplot(8, 10, (iconc) + (iconc-1)*8)
                timelfp = linspace(-0.5, 2, size(stats.Z{idx}, 2)); 
                plot(timelfp, mean(stats.Z{idx}+ii*yoffset, 1))
                hold on

                % Find inhalation 

                start_time = FVO(idx) - 0;
                end_time = FVO(idx) + 2;
                start_idx = searchsorted(SniffTime, start_time);
                end_idx = searchsorted(SniffTime, end_time);
                % Extract respiration traces
                
                postfvoinh = find(Prex>FVO(idx)); 
                inh1 = Prex(postfvoinh(1)); 
                inh2 = Prex(postfvoinh(2)); 

                respiration_trace_window = Sniff(start_idx:end_idx)+ ii*yoffset; 
                time_respiration = SniffTime(start_idx:end_idx)-inh1;
                %plot(time_respiration, respiration_trace_window, 'r')
                title([odorlabels{iodor}, ' conc = ', num2str(concentration) ])
            
                
                               
            end
            
            %imagesc(total_activity, [0, 1])
            % draw white line at x = 500
            hold on;
            plot( [0, 0], [0, size(total_activity, 1)],'k--', 'LineWidth', 2);
            % put xticks 
            %xticks_pos = linspace(1, size(lfp, 2), 5); 
            %set(gca, 'XTick', xticks_pos, 'XTickLabel', round(xticks_time, 2));
            %set(gca, 'YTick', yticks_pos, 'YTickLabel', round(yticks_freq, 2));
            xlabel('Time (s)');
            ylabel('Trial repetiton');
            
        
        end
    end
    clear stats
end

