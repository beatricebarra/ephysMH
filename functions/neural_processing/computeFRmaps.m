function [FRMAT , meanFRMAT, meanLATMAT,  LABELS, tvect] = computeFRmaps(spiketime_trains, goodunits, trials_numbers,  sig, mint, maxt)
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
FRMAT = cell(length(odorlabels),length(trials_numbers.(odorlabels{1}).concs)); 
LABELS.odor = cell(1, length(odorlabels)); % array (list of odor names)
LABELS.conc = cell(length(odorlabels), length(trials_numbers.(odorlabels{1}).concs));  % matrix (odors have different concentrations
iiu = 0; 
tvect = linspace(mint, maxt, (maxt-mint)/sig); 
for iodor = 1: length(odorlabels)
    for iconc = 1 : length(trials_numbers.(odorlabels{iodor}).concs)
        FRMAT{iodor, iconc} = zeros(length(goodunits), length(tvect)); 
        
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
end

for iodor = 1: length(odorlabels)
    meanFRMAT{iodor} = []; 
    meanLATMAT{iodor} = []; 
    for iconc = 1 : length(trials_numbers.(odorlabels{iodor}).concs)
        for iu = 1 : length(Rast{iodor, iconc})
            if isempty( Rast{iodor, iconc}{iu})
                FR = zeros(1, length(tvect)); 
            else
                if sum(cellfun(@(x) isempty(x), Rast{iodor, iconc}{iu})) == length(Rast{iodor, iconc}{iu})
                    FR = zeros(1, length(tvect)); 
                else
                    [FR,t, err] = compute_FR_from_raster(Rast{iodor, iconc}{iu} , sig, mint, maxt, tvect ); 
                end
            end
             % Attention this is the chronux findpeaks
            [meanFRMAT{iodor}(iu, iconc), meanLATMAT{iodor}(iu, iconc)] = max(FR); % for now with the max, but should probably do findpeaks or average in a window
            FRMAT{iodor, iconc}(iu, :)=  FR; 
        end
    end
end
%tvect = linspace(mint, maxt, (maxt-mint)/sig); 
%FR = zeros(1, length(tvect)); 
       
%tvect = linspace(mint, maxt, (maxt-mint)/sig); 
%[FR,t, err] = compute_FR_from_raster(Rast , sig, mint, maxt, tvect ); 
%FRMAT{iodor, iconc}(iiu, :)=  FR; 
                
    %end
%end



rmpath(genpath('/Users/barrab01/Documents/Repos/chronux_2_12'));
end