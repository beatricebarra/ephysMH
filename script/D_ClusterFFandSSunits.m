%% Computation of psths and relative ordering by peak latency

addpath(genpath('/Users/barrab01/Documents/Repos/matnwb'));
addpath(genpath('/Users/barrab01/Documents/Repos/ephysMH'));
addpath(genpath('/Users/barrab01/Documents/Repos/chronux_2_12'));
clc


clear all


% Mixing all mice 
filenames = {...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/UpdatedDatawWaveforms/040118__24_11_15.nwb',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/UpdatedDatawWaveforms/040118__24_11_18.nwb',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/UpdatedDatawWaveforms/040136__24_11_11.nwb',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/UpdatedDatawWaveforms/040136__24_11_12.nwb',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/UpdatedDatawWaveforms/040108__24_11_14.nwb', ...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/UpdatedDatawWaveforms/040108__24_11_16.nwb', ...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/UpdatedDatawWaveforms/040116__24_11_13.nwb', ...
    }

% savefilenames = {...
%     '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_11/Extracted_040136__24_11_11.mat',...
%     '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_12/Extracted_040136__24_11_12.mat',...
%     '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040116/24_11_13/Extracted_040116__24_11_13.mat',...
%     '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_14/Extracted_040108__24_11_14.mat',...
%     '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/12_11_15/Extracted_040118__24_11_15.mat',...
%     '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_16/Extracted_040108__24_11_16.mat',...
%     '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/24_11_18/Extracted_040118__24_11_18.mat',...
%     }
analogFS = 1000; 
dt = 1/analogFS; 
odorlabels = {'ethyltiglate', 'ethylbutyrate', 'acetophenone', 'heptanone'}; 
%%
ifilename  = 6

filename = filenames{ifilename}; 
nwbFile = nwbRead(filename);
%
spikes = nwbFile.units.spike_times.data(:); 
% Quick PCA on waveform features

waveform  = nwbFile.units.vectordata.get('waveform').data(:);
timevect = linspace(0, dt*length(waveform), length(waveform)); 

waveform_mat = reshape(waveform,  900, length(waveform)/900);

figure
plot( waveform, 'k')
hold on
plot(nwbFile.units.spike_times.data(:)/1000, ones(size(nwbFile.units.spike_times.data(:))), 'o')
% Waveform features
peak_to_valley  = nwbFile.units.vectordata.get('waveform_peak_to_valley').data(:);
waveform_asymetry  = nwbFile.units.vectordata.get('waveform_asymetry').data(:); 
waveform_half_width = nwbFile.units.vectordata.get('waveform_half_width').data(:); 
waveform_peak_trough_ratio = nwbFile.units.vectordata.get('waveform_peak_trough_ratio').data(:); 
waveform_recovery_slope = nwbFile.units.vectordata.get('waveform_recovery_slope').data(:); 
waveform_repolarization_slope = nwbFile.units.vectordata.get('waveform_repolarization_slope').data(:); 

data = [peak_to_valley, waveform_asymetry, waveform_half_width, waveform_peak_trough_ratio, waveform_recovery_slope, waveform_repolarization_slope]; 
rawwaveforms = waveform_mat'; 
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(rawwaveforms); 
figure
%plot3(SCORE(:,1), SCORE(:,2),  SCORE(:,3),'o')
cols = ['r', 'b', 'g', 'k']
[waveformID] = kmeans(SCORE(:,1:5), 3); 
for id = unique(waveformID)'
    idx = find(waveformID == id); 
    plot3(SCORE(idx,1), SCORE(idx,2),  SCORE(idx,3),'o', 'Color', cols(id))
    hold on
end

figure
for id = unique(waveformID)'
    idx = find(waveformID == id); 
    plot(-waveform_peak_trough_ratio(idx), waveform_recovery_slope(idx), 'o', 'Color', cols(id))
    hold on
end


figure
for id = unique(waveformID)'
    idx = find(waveformID == id); 
    meanwave = mean(rawwaveforms(idx, :)); 
    hold on
    plot(meanwave)
    
end

figure
hist(waveform_half_width*2, 1000)

idx_fs = find(waveform_half_width*2 <0.0004); 
idx_ss = intersect(find(waveform_half_width*2 >=0.0004), find(waveform_half_width*2 <=0.001)); 

figure

meanwavefs = mean(waveform_mat(:, idx_fs)'); 
meanwavess = mean(waveform_mat(:, idx_ss)'); 
hold on

stdshadesem(waveform_mat(:, idx_fs)', 0.2, 'r')
stdshadesem(waveform_mat(:, idx_ss)', 0.2, 'b')

    

%% Processs all mice together PCx 

sig = 0.010; 
mint = 0; 
maxt = 0.5; 
peak_latency = []; 
imouse = []; 
all_mice_condition_table= [];
all_mice_spikes_table = []; 
allmice_FRmat = []
for ifilename = 1: 1%length(filenames)
    
    filename = filenames{ifilename}; 
    nwbFile = nwbRead(filename);
    % Retrieve start and stop times
    [Odors, Concentration, FVO, FVC, Sniff, SniffTime, Prex, Postx, spikes, trials_start_stop_time] = extract_variables(nwbFile); 
    % Process sniff 
    all_odors = unique(Odors); 
    Sniffdt = mean(diff(SniffTime)); 
    SniffFS = 1/Sniffdt; 
    fSniff = lowPass(Sniff, round(1/Sniffdt), 10, 3); 
    baselineSniff = lowPass(Sniff, round(1/Sniffdt),0.8, 2); 
    % Quality control
    [OB_units_idx, PC_units_idx] = units_quality_check(nwbFile); 
    % Extract spikes
    % spikes table has #units rows and #conditions columns
    [spikes_table, count_table, condition_table] = extract_spikes(nwbFile, spikes, Odors, Concentration, Prex, Postx, FVO,  PC_units_idx, all_odors); 
    all_mice_condition_table = [all_mice_condition_table; condition_table]; 
    all_mice_spikes_table = [all_mice_spikes_table; spikes_table]; 
    % Extract peak latency of PSTH
    [FRmatthismouse, peak_latencythismouse] = extract_peak_latency(spikes_table, sig, mint, maxt); 
    peak_latency= [peak_latency; peak_latencythismouse]; 
    % Initialize only if it's the first mouse...
    if ifilename ==1
        allmice_FRmat = cell(1, size(peak_latency, 2)); 
    end
    % .. and then fill matrix at every mouse
    for icondition = 1 : size(peak_latency, 2)
        allmice_FRmat{icondition}= [allmice_FRmat{icondition}; FRmatthismouse{icondition}]; 
    end

    disp([ 'Dimensions: spikes_table : ', num2str(size(spikes_table, 1)), ', latency: ', num2str(size(peak_latencythismouse, 1)), ', FR: ', num2str(size(FRmatthismouse{1}, 1))])
    mouse_label = ifilename*ones(size(peak_latency, 1), 1); 
    imouse = [imouse; ifilename*ones(size(peak_latencythismouse, 1), 1)]; 
    %save(savefilenames{ifilename}, 'peak_latencythismouse', 'FRmatthismouse', 'spikes_table', 'condition_table', 'mouse_label', "-v7.3")
end
%savefilenameallmice = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/all_mice_dataset.mat'; 
%save(savefilenameallmice, 'peak_latency', 'allmice_FRmat', 'all_mice_spikes_table', 'all_mice_condition_table', 'imouse', "-v7.3")

