addpath(genpath('/Users/barrab01/Documents/Repos/matnwb'));
addpath(genpath('/Users/barrab01/Documents/Repos/ephysMH'));
addpath(genpath('/Users/barrab01/Documents/Repos/chronux_2_12'));

folder = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_14/'
filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_14/040108__24_11_14.nwb'

%folder = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_16/'
%filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_16/040108__24_11_16.nwb'

%folder = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040116/24_11_13/'
%filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040116/24_11_13/040116__24_11_13.nwb'

%folder = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/12_11_15/'
%filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/12_11_15/040118__24_11_15.nwb'

%folder = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/24_11_18/'
%filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/24_11_18/040118__24_11_18.nwb'

folder = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_11/'
filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_11/040136__24_11_11.nwb'

savefilename = fullfile(folder, '/reduced_dataset.mat')
savedistrialname = fullfile(folder, 'discarded_trials.mat')
savetrialnumbers = fullfile(folder, 'trials_numbers.mat')

nwbFile = nwbRead(filename);
analogFS = 1000; 
dt = 1/analogFS; 
odorlabels = {'ethyltiglate', 'ethylbutyrate', 'acetophenone', 'heptanone'}; 

%% Small initial script 
figure
% Retrieve start and stop times
trials_start_stop_time = [nwbFile.intervals.get('odor_trials').start_time.data(:), ...
                        nwbFile.intervals.get('odor_trials').stop_time.data(:)]

Odors = nwbFile.intervals.get('odor_trials').vectordata.get('pins').data(:); 
Concentration = nwbFile.intervals.get('odor_trials').vectordata.get('concentration').data(:); 
% Sniff data
Sniff = nwbFile.acquisition.get('respiration_data').data(:);  
SniffTime = nwbFile.acquisition.get('respiration_data').timestamps(:);  
Sniffdt = mean(diff(nwbFile.acquisition.get('respiration_data').timestamps(:))); 
fSniff = lowPass(Sniff, round(1/Sniffdt), 10, 3); 
baselineSniff = lowPass(Sniff, round(1/Sniffdt),0.8, 2); 
SniffFS = 1/mean(diff(nwbFile.acquisition.get('respiration_data').timestamps(:))); 



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
save(savetrialnumbers,"trials_numbers","-v7.3")


%%
% Params of experiment
plot_time = 2; 
Nodors = length(unique(Odors)); 
Nconc = 8; 
trial_found_table = zeros(Nodors, Nconc); 
for itrial = 1 : size(trials_start_stop_time, 1)
    itrial
    % Create raster of this trial
    %subplot(5, 5, itrial)
    spiketime_trains{itrial} = struct(); 
    % Odor 
    all_odors = unique(Odors);
    thisodor = Odors(itrial); 
    thisodor_idx = find(all_odors == thisodor); 
    idx_odor = find(Odors ==thisodor); 
    
    % .. and concentration 
    all_conc = unique(Concentration(idx_odor)); 
    thisconc = Concentration(itrial); 
    thisconc_idx = find(all_conc == thisconc); 
    spiketime_trains{itrial}.odor =  [thisodor, thisodor_idx]; 
    spiketime_trains{itrial}.conc  = [thisconc, thisconc_idx];
    %if trial_found_table(thisodor_idx, thisconc_idx) == 0 
    
    % Flag that you already plotted
    %trial_found_table(thisodor_idx, thisconc_idx) = 1; 
    %subplot(Nodors, Nconc, (thisodor_idx-1)*Nconc + thisconc_idx)

    % Find sniff snip (filtered and baseline too)
    [minv, idxSniff] = min(abs(SniffTime - trials_start_stop_time(itrial, 1))); 
    SniffSnip = Sniff(idxSniff : idxSniff + floor(plot_time/ Sniffdt)); 
    fSniffSnip = fSniff(idxSniff : idxSniff + floor(plot_time/ Sniffdt));
    basSniffSnip = baselineSniff(idxSniff : idxSniff + floor(plot_time/ Sniffdt));
    SniffTimeSnip = SniffTime(idxSniff : idxSniff + floor(plot_time/ Sniffdt)); 
    spiketime_trains{itrial}.Sniff = SniffSnip; 
    spiketime_trains{itrial}.fSniff = fSniffSnip; 
    spiketime_trains{itrial}.bSniff = basSniffSnip; 
    spiketime_trains{itrial}.SniffTime = SniffTimeSnip; 
    
    

    %plot(SniffTimeSnip, (SniffSnip*100)+900, 'k')
    %hold on
    %title([num2str(Odors(itrial)), '_', num2str(Concentration(itrial))])

    spiketime_trains{itrial}.units = {}; 
    previous_unit_idx = 1; 
    for iunit = 1 : max(nwbFile.units.vectordata.get('cluster_id').data(:))-1 
        thisunit_idx =  previous_unit_idx; 
        nextunit_idx =  nwbFile.units.spike_times_index.data(iunit);  
        if strcmp(nwbFile.units.vectordata.get('location').data{iunit}{1}(1:2), 'OB')
            col = 'b';
            type = '1'; 
        elseif strcmp(nwbFile.units.vectordata.get('location').data{iunit}{1}(1:2), 'PC')
            col = 'r';
            type = '2'; 
        end
        % Filter for quality 
        if strcmp(nwbFile.units.vectordata.get('SI_quality').data{iunit}{1}, 'good') && strcmp(nwbFile.units.vectordata.get('quality').data{iunit}{1}, 'good')
            allspikes_unit = nwbFile.units.spike_times.data(thisunit_idx: nextunit_idx); 
            trial_spikes_unit_idx = intersect(find(allspikes_unit>trials_start_stop_time(itrial, 1)), find(allspikes_unit<trials_start_stop_time(itrial, 1)+plot_time));  
            trial_spikes_times = allspikes_unit(trial_spikes_unit_idx); 
            spiketime_trains{itrial}.units{iunit} = trial_spikes_times;
            spiketime_trains{itrial}.unittype(iunit) = type;
            mv = ones(size(trial_spikes_times)); 
%             for spiketime = trial_spikes_times'
%                 plot([spiketime, spiketime], [double(iunit), double(iunit-0.8)], col)
%                 hold on
%             end
        end
        previous_unit_idx = nextunit_idx+1; 
    end
end

%% Find and save inhalation

discarded_trials = []; 
for itrial = 1 : length(spiketime_trains)
    
    % Filter Sniff 

    filtered_sniff = spiketime_trains{itrial}.fSniff;  
    hold on
   
    th = mean(spiketime_trains{itrial}.bSniff) - 1*std(baselineSniff); 
    odortime = spiketime_trains{itrial}.SniffTime(1);
    logicSniff = (filtered_sniff < th);
    try
        if logicSniff(1) == 0
            all_ones = find(logicSniff ==1); 
            spiketime_trains{itrial}.inhalation_idx = all_ones(1); 
            spiketime_trains{itrial}.inhalation_time = spiketime_trains{itrial}.SniffTime(all_ones(1)); 
        else
            gaps = diff(find(logicSniff)); 
            idxg = find(gaps>20); 
            idx = idxg(1)+gaps(idxg(1)); 
            spiketime_trains{itrial}.inhalation_idx = idx; 
            spiketime_trains{itrial}.inhalation_time = spiketime_trains{itrial}.SniffTime(idx); 
        end
        
        iter_idx = spiketime_trains{itrial}.inhalation_idx-1; 
        baseline_not_crossed = 1; 
        while baseline_not_crossed
            if filtered_sniff(iter_idx)> mean(spiketime_trains{itrial}.bSniff)
                baseline_not_crossed = 0; 
            else
                iter_idx = iter_idx-1; 
            end
        end
        spiketime_trains{itrial}.inhalation_idx = iter_idx; 
        spiketime_trains{itrial}.inhalation_time = spiketime_trains{itrial}.SniffTime(iter_idx);        
    catch
        spiketime_trains{itrial}.inhalation_idx = NaN; 
        spiketime_trains{itrial}.inhalation_time = NaN; 
        
        discarded_trials = [discarded_trials; itrial]; 
        
    end  
end




%% Save

save(savefilename,"spiketime_trains","-v7.3")
save(savedistrialname,"discarded_trials","-v7.3")
save(savetrialnumbers,"trials_numbers","-v7.3")


