function [PCIDX, exPCCL] = findExcitatoryUnits(spiketime_trains, trials_numbers, sig, mint, maxt, unittype)
    % Find excitatory units PCx through PCA
    % Plot PCA and kmeans analysis for PC and OB units
    addpath(genpath('/Users/barrab01/Documents/Repos/chronux_2_12'));

    odorlabels = fieldnames(trials_numbers);
    % Extract good PC and OB units
    goodPCunits = strfind(spiketime_trains{1}.unittype, unittype);
    % PC analysis first
    PCMAT = []; 
    for iunit = goodPCunits
        PCline = []; 
        labelmat = []; 
        for iodor = 1: length(odorlabels)
            for iconc = 1 : length(trials_numbers.(odorlabels{iodor}).concs)
                labelmat = [labelmat, [iodor; iconc]]; 
                try
                    Rast = {}; 
                    for ii = 1 : length(trials_numbers.(odorlabels{iodor}).trial_idxs{iconc})
                        itrial = trials_numbers.(odorlabels{iodor}).trial_idxs{iconc}(ii); 
                        idx = intersect(find(spiketime_trains{itrial}.units{iunit}> spiketime_trains{itrial}.inhalation_time-mint), ...
                                            find(spiketime_trains{itrial}.units{iunit}< spiketime_trains{itrial}.inhalation_time + maxt)); 
                        spikes = spiketime_trains{itrial}.units{iunit}(idx)-spiketime_trains{itrial}.inhalation_time; 
                        Rast{ii} = spikes;
                    end
                    tvect = linspace(mint, maxt, (maxt-mint)/sig); 
                    [FR,t, err] = compute_FR_from_raster(Rast , sig, mint, maxt, tvect ); 
                catch
                    FR = 0; 
                end
                meanFR = mean(FR); 
                PCline= [PCline,meanFR ]; 
            end
        end
        PCMAT= [PCMAT; PCline]; 
    end
    PCMATraw = PCMAT; 
    PCMAT = (zscore(PCMATraw')'); 
    clustercols = {'r', 'g', 'b', 'y'}; 
    [clusterPCCells, PCIDX] = cluster_units(PCMAT, clustercols); 
    for i = unique(PCIDX)'
        PCClSize(i) = length(find(PCIDX ==i)); 
    end
    [U, exPCCL] = max(PCClSize); 
    rmpath(genpath('/Users/barrab01/Documents/Repos/chronux_2_12'));

end