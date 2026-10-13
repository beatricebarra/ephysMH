
%  Bea's paths
addpath(genpath('/Users/barrab01/Documents/Repos/PiriformAnalysis/Dependencies'))
addpath(genpath('/Users/barrab01/Documents/Repos/chronux_2_12'))
addpath(genpath('/Users/barrab01/Documents/Repos/fieldtrip'))
addpath(genpath('/Users/barrab01/Documents/Repos/matnwb'))
rmpath('/Users/barrab01/Documents/Repos/fieldtrip/external/signal/')

filename_units = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_14/040108__24_11_14.nwb'
LFPfolder = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_14/LFP'; 

% Create save folder path
savefolder = [LFPfolder, '/processed/']; 

% Load LPF data
map_file_path = fullfile(LFPfolder, 'PC_LFP_channelMap.csv'); 
stamps_file_path = fullfile(LFPfolder, 'PC_LFP_timestamps.csv'); 
data_file_path = fullfile(LFPfolder, 'PC_lfp.dat'); 
% Read CSV map file
lfp_map = readmatrix(map_file_path);   % or readtable if you need table format

% Load stamps file (comma-delimited)
stamps = readmatrix(stamps_file_path); 

% Get dimensions
nSamples  = size(stamps, 1);
nChannels = numel(lfp_map);   % equivalent to flatten().shape[0]

% Memory-map the binary data file
m = memmapfile(data_file_path, ...
    'Format', {'int16', [nSamples, nChannels], 'lfp'}, ...
    'Writable', false);

% Access mapped data as:
lfp = m.Data.lfp;

% Read start and stop of trials
nwbFile = nwbRead(filename_units);
% Read important variables
trials_start_stop_time = [nwbFile.intervals.get('odor_trials').start_time.data(:), ...
                        nwbFile.intervals.get('odor_trials').stop_time.data(:)]
Odors = nwbFile.intervals.get('odor_trials').vectordata.get('pins').data(:); 
Concentration = nwbFile.intervals.get('odor_trials').vectordata.get('concentration').data(:); 
FVO = nwbFile.intervals.get('odor_trials').vectordata.get('fvo').data(:); 
FVC = nwbFile.intervals.get('odor_trials').vectordata.get('fvc').data(:); 
% Read Sniff features
Prex = nwbFile.analysis.get('respiratory_features').vectordata.get('prex').data(:); 
Postx = nwbFile.analysis.get('respiratory_features').vectordata.get('postx').data(:); 
Sniff = nwbFile.acquisition.get('respiration_data').data(:);  
SniffTime = nwbFile.acquisition.get('respiration_data').timestamps(:);  
Sniffdt = mean(diff(nwbFile.acquisition.get('respiration_data').timestamps(:))); 
fSniff = lowPass(Sniff, round(1/Sniffdt), 10, 3); 
baselineSniff = lowPass(Sniff, round(1/Sniffdt),0.8, 2); 
SniffFS = 1/mean(diff(nwbFile.acquisition.get('respiration_data').timestamps(:))); 

% Find trial times
events = []
for itrial = 1 : length(Odors)
    start_time = FVO(itrial) - 10;
    end_time = FVO(itrial) + 6;
    start_idx = searchsorted(SniffTime, start_time);
    end_idx = searchsorted(SniffTime, end_time);
    % Extract respiration traces
    postfvoinh = find(Prex>FVO(itrial)); 
    inh1 = Prex(postfvoinh(1)); 
    inh2 = Prex(postfvoinh(2)); 
    events = [events;inh1]
    
end



% 
% figure
% hold on
% 
% for itrial = 1 : size(trials_start_stop_time, 1)-200
%     start_trial = trials_start_stop_time(itrial,1);
%     stop_trial = trials_start_stop_time(itrial,2); 
% 
%     [~, idx_start] = min(abs(stamps - start_trial)); 
%     [~, idx_stop] = min(abs(stamps - stop_trial)); 
%     disp([idx_start, idx_stop])
% 
%     LFP_this_trial = lfp(idx_start:idx_stop, 1)';
%     plot(LFP_this_trial)
%     pause
% 
% end

%% Set the parameters for spectrum computation

para = [];
% previous parameters
para.freq = logspace(log10(0.5),log10(200),102); % center frequencies for filtering
para.bandwidth = logspace(log10(1),log10(30),102); % bandwidth of each center frequency
% trying this with lower hi-freq cutoff
% log scale
%     para.freq = logspace(log10(0.5),log10(100),51); % center frequencies for filtering
%     para.bandwidth = logspace(log10(1),log10(15),51); % bandwidth of each center frequency
% log scale - higher freq
%     para.freq = logspace(log10(0.5),log10(200),102);
%     para.bandwidth = logspace(log10(1),log10(30),102);
% lin scale
%     para.freq = linspace(0.5,200,150);
%     para.bandwidth = linspace(0.5,20,150);
% lin scale - higher freq
% % Log Freq Low Freq
%     para.freq = logspace(log10(0.5),log10(30),70); % center frequencies for filtering
%     para.bandwidth = logspace(log10(1),log10(10),70); % bandwidth of each center frequency
% Low Freq Linear
%     para.freq = (0.5:0.5:29);
%     para.bandwidth = ones(1,length(para.freq));
%     % High Freq Linear
%     para.freq = (30:2.5:200);
%     para.bandwidth = ones(1,length(para.freq)) * 5;

para.srate = 1000; % sampling rate, in Hz
idx_event = []
for event = events'
    [~, idx] = min(abs(stamps - event)); 
    idx_event = [idx_event; idx]; 
end

para.events = num2cell(idx_event); % in samples
para.epoch = [-0.5, 2]; % in seconds, relative to events
para.baseline = [-0.3, 0]; % in seconds, relative to events
% para.baseline = [-1, 0]; % in seconds, relative to events
para.smooth_win = 0.01; % temporal smoothing, in secondsj
para.nbpermuts = 10000; % number of permutations,
para.measure = 'amplitude'; % or 'power'
para.baseline_condition = []; % which condition to use as baseline (empty for just before each of multiple conditions)
para.trl = 'mean'; % 'mean' or 'median' value across trials

%% Compute Spectrogram for for each channel
h = waitbar(0, 'Computation of Spectrogram...');

for ich = 1 : size(lfp, 2)
    % Display wait bar
    waitbar(ich/size(lfp, 2), h, sprintf('Progress: %d%%', ich/size(lfp, 2)));
    
    % Compoute spectrogram amplitude on the whole downsampled signal
    [amp, baseline_condition] = ComputeAmplitude(lfp(:,ich), para); 
    % Compute baseline correction
    [resp, resp_percentchange, nbevs] = correctBaseline(amp, baseline_condition, para); 
    % Perform permutation and compute Z 
    [stats] =performPermutation(amp, resp, nbevs, para, lfp(:,ich));
    % Save file
    pathsave = fullfile(savefolder, ['stats_ch', num2str(ich), '.mat']);
    save(pathsave, 'stats', '-v7.3');
    %[stats, epoch_time] = PermutAmpChange_Commonbaseline(lfp(:,1), para); % this function is in   /Users/barrab01/Documents/Repos/PiriformAnalysis/Dependencies
end


%%
% for ich = 1 : size(lfp, 2)
%     pathsave = fullfile(savefolder, ['stats_ch', num2str(ich), '.mat'])
% end
%% Plot 

% figure
% nOdors = length(unique(Odors)); 
% nConcentrations = 8; 
% iodor = 0; 
% for odor = unique(Odors)'
%     iodor = iodor +1; 
%     idx_odor =find( Odors == odor); 
%     iconc = 0; 
%     for conc = unique(Concentration(idx_odor))'
%         iconc = iconc + 1; 
%         idx_conc = find( Concentration == conc); 
%         % intersect idx 
%         idx_trials = intersect(idx_odor, idx_conc); 
%         % SUm spetrum of different trials
%         AvgSpectrum = zeros(size(stats.Z{1})); 
% 
%         for idx_t = idx_trials'
%             AvgSpectrum = AvgSpectrum + stats.Z{idx_t}; 
%         end
%         AvgSpectrum = AvgSpectrum./length(idx_trials); 
%         subplot(nOdors, nConcentrations, ((iodor-1)*nConcentrations + iconc))
%         imagesc(AvgSpectrum)
%     end
% end