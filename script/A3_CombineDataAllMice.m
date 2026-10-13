%% Computation of psths and relative ordering by peak latency

addpath(genpath('/Users/barrab01/Documents/Repos/matnwb'));
addpath(genpath('/Users/barrab01/Documents/Repos/ephysMH'));
addpath(genpath('/Users/barrab01/Documents/Repos/chronux_2_12'));
clc

clear all

folder = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_14/';
filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_14/040108__24_11_14.nwb';

folder ='/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_16/';
filename ='/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_16/040108__24_11_16.nwb';


% Mixing all mice 
filenames = {...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_11/040136__24_11_11.nwb',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_12/040136__24_11_12.nwb',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040116/24_11_13/040116__24_11_13.nwb',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_14/040108__24_11_14.nwb',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/12_11_15/040118__24_11_15.nwb',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_16/040108__24_11_16.nwb',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/24_11_18/040118__24_11_18.nwb',...
    }

savefilenames = {...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_11/Extracted_040136__24_11_11.mat',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_12/Extracted_040136__24_11_12.mat',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040116/24_11_13/Extracted_040116__24_11_13.mat',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_14/Extracted_040108__24_11_14.mat',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/12_11_15/Extracted_040118__24_11_15.mat',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_16/Extracted_040108__24_11_16.mat',...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/24_11_18/Extracted_040118__24_11_18.mat',...
    }
analogFS = 1000; 
dt = 1/analogFS; 
odorlabels = {'ethyltiglate', 'ethylbutyrate', 'acetophenone', 'heptanone'}; 

%% Processs all mice together PCx 

sig = 0.010; 
mint = 0; 
maxt = 0.5; 
peak_latency = []; 
imouse = []; 
all_mice_condition_table= [];
all_mice_spikes_table = []; 
allmice_FRmat = []
for ifilename = 1: length(filenames)
    
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
    save(savefilenames{ifilename}, 'peak_latencythismouse', 'FRmatthismouse', 'spikes_table', 'condition_table', 'mouse_label', "-v7.3")
end
savefilenameallmice = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/all_mice_dataset.mat'; 
save(savefilenameallmice, 'peak_latency', 'allmice_FRmat', 'all_mice_spikes_table', 'all_mice_condition_table', 'imouse', "-v7.3")

