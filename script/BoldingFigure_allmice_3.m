% Make Kevin Bolding's figure for each session
clear all
close all

session_paths = {
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_14', ...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040108/24_11_16', ...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040116/24_11_13', ...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/12_11_15', ...
    '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/24_11_18', ...
    %'/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_11', ...
    %'/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_12', ...
    }; 
odorlabels = {'ethyltiglate', 'ethylbutyrate', 'acetophenone', 'heptanone'}; 

%% Loop over all sessions and load
sessions = {}
for ipath = 1 : length(session_paths)
    sessions{ipath} = load(fullfile(session_paths{ipath}, 'reduced_dataset.mat')); 
    temp =  load(fullfile(session_paths{ipath}, 'trials_numbers.mat')); 
    sessions{ipath}.trial_numbers = temp.trials_numbers; 
    
end
 


imouse = 3; 
%% Raster plots
figure


goodPCunits = strfind(sessions{imouse}.spiketime_trains{1}.unittype, '2');
goodOBunits = strfind(sessions{imouse}.spiketime_trains{1}.unittype, '1');
sig = 0.020; 
mint = -0.3; 
maxt = 0.3; 
goodunits =goodPCunits; 

[Rast] = computeRasters(sessions{imouse}.spiketime_trains, goodunits, sessions{imouse}.trial_numbers, sig, mint, maxt); % Extract Rasters
mycolors = hot(15); 
plotRasters_withrepetitions(Rast, odor_idx, conc_idx, mycolors)
%%
% Params
i = 3; 
sig = 0.010; 
mint = -0.1; 
maxt = 0.2; 

goodPCunits = strfind(sessions{i}.spiketime_trains{1}.unittype, '2');
goodOBunits = strfind(sessions{i}.spiketime_trains{1}.unittype, '1');

colormap("pink")
goodunits =goodPCunits; 

goodunits = goodOBunits; 
[FRMAT , LABELS, tvect] = computeFRmaps(sessions{i}.spiketime_trains, goodunits, sessions{i}.trial_numbers, sig, mint, maxt); % Extract FR maps

figure
Nstd = 1.5; 
for i = 1:4
    for j = 1 : 8
        ii = (i-1)*8 + j; 
        
        subplot(4, 8, ii)
        
        [sortedFRMAT,Isort, sortedv]  = sort_cells_by_peak(FRMAT{i, j}, tvect, Nstd); 
        nonsortedmat = FRMAT{i, j}; 
        
        imagesc(sortedFRMAT, [0, 90])
        hold on
        idx_inh = floor(size(sortedFRMAT, 2)*abs(-mint/(-mint + maxt))); 
        plot([idx_inh, idx_inh], [0, size(sortedFRMAT, 1)], '--w', 'Linewidth', 2) 
        %[M, I] = max(sortedFRMAT, [], 2); 
        p = plot(sortedv, flipud(linspace( 1, size(sortedFRMAT, 1), size(sortedFRMAT, 1))), 'ow', 'MarkerFaceColor', 'w', 'MarkerSize', 2)
        hold off
        alpha(p, 0.1)
        
       
        title([LABELS.odor{i}, '_',num2str(LABELS.conc{i,j})])
        colormap("pink")
        if j == 6
            subplot(4, 8, ii-4)
            hold on
            p = plot(sortedv, flipud(linspace( 1, size(sortedFRMAT, 1), size(sortedFRMAT, 1))), 'oy', 'MarkerFaceColor', 'y', 'MarkerSize', 2)
        end

        
    end
end



%% 
for isession = 1 : length(sessions)
    % For each session: 
    spiketime_trains = sessions{isession}.spiketime_trains; 
    % Find OB and PC units (from first trial): list should be the same for
    % all trial
    goodPCunits = strfind(spiketime_trains{1}.unittype, '2');
    goodOBunits = strfind(spiketime_trains{1}.unittype, '1');



end
%% 

% I have to write these functions : 
% - Create FR map 
% - Create histograms (population) with gaussian components
% - divie cells in inhibitory and excitatory (you have to understand how bolding does it) 



