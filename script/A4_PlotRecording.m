% Check recording of a specific mouse 

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

analogFS = 1000; 
dt = 1/analogFS; 
odorlabels = {'ethyltiglate', 'ethylbutyrate', 'acetophenone', 'heptanone'}; 
intensity_file  = '/Users/barrab01/Documents/Repos/ephysMH/perceived_intensity_curves.xlsx'; 

%% Processing parameters

sig = 0.010; % Window of gaussian kernel
mint = 0; % from inhalation start, time at which plot window starts
maxt = 0.5; % from inhalation start, time at which plot window ends

%% 
allI = [];
allC = []; 
allISI = [];
for ifilename = 1:length(filenames)
    
    % Read the file 
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
    %[spikes_table, count_table, condition_table] = extract_spikes(nwbFile, spikes, Odors, Concentration, Prex, Postx, FVO,  PC_units_idx, all_odors);

    % Visualize respitation traces and compute inter-sniff- time 
    ISI = plot_session_respiration(Odors, Concentration, FVO, FVC, fSniff, SniffTime,Prex)
    
    % Work on sniff 
    figISIMouse = figure
    plot_session_SniffIntervals(figISIMouse, Odors, Concentration, ISI);
    
    io = 0; 
    for o = unique(Odors)'
        conc_this_odor = unique(Concentration(find(Odors == o))); 
        allconc = unique(conc_this_odor); 
        io = io+1; 
        ic = 1; 
        for c = allconc(2:end)'
            ic = ic +1; 
             % Compute intensity of each concentration 
            iodor = o-4; 
            [I, pC ]= computeIntensity(intensity_file, [double(iodor), c]); 
            realI(io, ic, :) = [ones(1, length(squeeze(ISI(io, ic, :))))*I]; 
            realC(io, ic, :) = [ones(1, length(squeeze(ISI(io, ic, :))))*pC]; 
        end
    end
    allI = [allI; realI(:)]; 
    allC = [allC; realC(:)];
    allISI = [allISI; ISI(:)]; 
end


figure
% Compute regression intensity - sniff 
% Ensure column vectors
my_linear_regression( allI, allISI)

% Compute regression concentration (log10(ppm)) - sniff 
my_linear_regression( allC, allISI)


%% Random forest for sniff vs Intensity 

my_random_forest(allI, allISI)

my_random_forest(allC, allISI)

