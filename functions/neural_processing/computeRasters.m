function [Rast] = computeRasters(spiketime_trains, goodunits, trials_numbers,  sig, mint, maxt)
% Inputs: 
% -spiketime_trains : cell array containing sniff and unit spikes for each trial
% -goodunits: selected units
% -trials_numbers: suddivision of trials per condition
% Outputs: 
% -FRMAT : cell array containing matrixes or firing rates for each unit (not sorted by activation time)
% -LABELS: labels for condition corresponding to FRMAT

addpath(genpath('/Users/barrab01/Documents/Repos/chronux_2_12'));


% PC analysis first
odorlabels = fieldnames(trials_numbers); 
LABELS.odor = cell(1, length(odorlabels)); % array (list of odor names)
LABELS.conc = cell(length(odorlabels), length(trials_numbers.(odorlabels{1}).concs));  % matrix (odors have different concentrations
for iodor = 1: length(odorlabels)
    for iconc = 1 : length(trials_numbers.(odorlabels{iodor}).concs)
        repcount{iodor, iconc} = 0; 
        Rast{iodor, iconc}= cell(1, length(goodunits)); 
        for iu = 1 : length(goodunits)
            Rast{iodor, iconc}{iu} = {}; 
        end
    end
end

for itrial = 1 : length(spiketime_trains)
    iodor = spiketime_trains{itrial}.odor(2); 
    iconc = spiketime_trains{itrial}.conc(2); 
    repcount{iodor, iconc} = repcount{iodor, iconc} +1; 
    LABELS.odor{iodor} = odorlabels{iodor};      
    LABELS.conc{iodor, iconc} = trials_numbers.(odorlabels{iodor}).concs(iconc); 
    iiu = 0; 
    for iunit = goodunits
        iiu = iiu + 1; 
        
        idx = intersect(find(spiketime_trains{itrial}.units{iunit}> spiketime_trains{itrial}.inhalation_time+mint), ...
                            find(spiketime_trains{itrial}.units{iunit}< spiketime_trains{itrial}.inhalation_time + maxt)); 
        spikes = spiketime_trains{itrial}.units{iunit}(idx)-spiketime_trains{itrial}.inhalation_time; 
        Rast{iodor, iconc}{iiu}{repcount{iodor, iconc}} = spikes;
            
            
    end
rmpath(genpath('/Users/barrab01/Documents/Repos/chronux_2_12'));

end