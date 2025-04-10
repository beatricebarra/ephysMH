function [Odors, Concentration, FVO, FVC, Sniff, SniffTime, Prex, Postx, spikes, trials_start_stop_time] = extract_variables(nwbFile)
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
    % Sniff features
    Prex = nwbFile.analysis.get('respiratory_features').vectordata.get('prex').data(:); 
    Postx = nwbFile.analysis.get('respiratory_features').vectordata.get('postx').data(:); 
    % Spikes
    spikes = nwbFile.units.spike_times.data(:); 
end