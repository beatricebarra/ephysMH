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
sig = 0.010; % 10 instead of 20 look at methods section
mint = -0.2; 
maxt = 0.6; 
Nstd = 1.5; 
%% Loop over all sessions and load
sessions = {}
for ipath = 1 : length(session_paths)
    sessions{ipath} = load(fullfile(session_paths{ipath}, 'reduced_dataset.mat')); 
    temp =  load(fullfile(session_paths{ipath}, 'trials_numbers.mat')); 
    sessions{ipath}.trial_numbers = temp.trials_numbers; 
    
end

%% Quality filtering for units 


%% 
figure
for imouse = 1 : 5

    
    goodPCunits = strfind(sessions{imouse}.spiketime_trains{2}.unittype, '2');
    goodOBunits = strfind(sessions{imouse}.spiketime_trains{2}.unittype, '1');
    %[PCIDX, exPCCL] = findExcitatoryUnits(sessions{imouse}.spiketime_trains, sessions{imouse}.trial_numbers, sig, mint, maxt, '2'); % PC = 2
    %idx_ex = find(PCIDX ==exPCCL); 
    %goodunits =goodPCunits(idx_ex); 
    goodunits  = goodOBunits; 
    
    % Compute firing rate maps
    [FRMAT , meanFRmat, meanLATMAT, LABELS, tvect] = computeFRmaps(sessions{imouse}.spiketime_trains, goodunits, sessions{imouse}.trial_numbers, sig, mint, maxt); % Extract FR maps
    TEMP = []; 
    subplot(1,5,(imouse) )
    hold on
    for iodor = 1 : length(meanFRmat)  
        TEMP = [TEMP, meanFRmat{iodor}]; 
    end
    
    imagesc(TEMP)
    hold on
    for iodor = 1 : length(meanFRmat)
        plot([(iodor-1)*8+0.5, (iodor-1)*8+0.5], [0, size(TEMP, 1)], 'r', 'Linewidth', 0.8)
        hold on
        myxticks = [myxticks, (iodor-1)*8+1]; 
    end
    xticks(myxticks)
    xticklabels(odorlabels)
    % Compute raster plot
    %[Rast] = computeRasters(sessions{imouse}.spiketime_trains, goodunits, sessions{imouse}.trial_numbers, sig, mint, maxt); % Extract Rasters
    % Compute gaussian 
    %[MeanGaussians, x, UpperboundGaussians, LowerboundGaussians] = computeMeanGaussians(sessions{imouse}.spiketime_trains, sessions{imouse}.trial_numbers, goodunits, 0, 0.6, 3)
end


%% Produce a figure for each mouse
for imouse = 1 : 5

    
    goodPCunits = strfind(sessions{imouse}.spiketime_trains{2}.unittype, '2');
    goodOBunits = strfind(sessions{imouse}.spiketime_trains{2}.unittype, '1');
    %[PCIDX, exPCCL] = findExcitatoryUnits(sessions{imouse}.spiketime_trains, sessions{imouse}.trial_numbers, sig, mint, maxt, '2'); % PC = 2
    %idx_ex = find(PCIDX ==exPCCL); 
    %goodunits =goodPCunits(idx_ex); 
    goodunits  = goodPCunits; 
    
    % Compute firing rate maps
    [FRMAT , meanFRmat, LABELS, tvect] = computeFRmaps(sessions{imouse}.spiketime_trains, goodunits, sessions{imouse}.trial_numbers, sig, mint, maxt); % Extract FR maps
    
    % Compute raster plot
    [Rast] = computeRasters(sessions{imouse}.spiketime_trains, goodunits, sessions{imouse}.trial_numbers, sig, mint, maxt); % Extract Rasters
    % Compute gaussian 
    [MeanGaussians, x, UpperboundGaussians, LowerboundGaussians] = computeMeanGaussians(sessions{imouse}.spiketime_trains, sessions{imouse}.trial_numbers, goodunits, 0, 0.6, 3)
    
    % Plot everything in the same figure
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    concentration_idx = [1, 2, 6]; 
    mycolors = hot(12); 

    concentration_colors = parula(8); 
    for row = 1:3
        for col = 1:3
            % Function 1 (columns 1-3)
            subplot(3,9, (row-1)*9 + col);
            plotRaster_withrepetitions( Rast{row, concentration_idx(col)}, mycolors)
            title(sprintf('Odor %d, Conc %d', row, concentration_idx(col)));
        
            % Function 2 (columns 4-6)
            subplot(3,9, (row-1)*9 + col + 3);
            [sortedFRMAT,Isort, sortedv]  = sort_cells_by_peak(FRMAT{row, concentration_idx(col)}, tvect, Nstd); 
            imagesc(sortedFRMAT, [0, 50])
            hold on
            idx_inh = floor(size(sortedFRMAT, 2)*abs(-mint/(-mint + maxt))); 
            plot([idx_inh, idx_inh], [0, size(sortedFRMAT, 1)], '--w', 'Linewidth', 2) 
            %[M, I] = max(sortedFRMAT, [], 2); 
            p = plot(sortedv, flipud(linspace( 1, size(sortedFRMAT, 1), size(sortedFRMAT, 1))), 'ow', 'MarkerFaceColor', 'w', 'MarkerSize', 2)
            hold off
            %set(gca, 'xticks', linspace(1, size(sortedFRMAT, 2), 5),'xtickslabel', linspace(mint, maxt, 5))
            xticks(linspace(1, size(sortedFRMAT, 2), 5))
            xticklabels(linspace(mint, maxt, 5))
            alpha(p, 0.1)
            title(sprintf('Odor %d, Conc %d', row, concentration_idx(col)));
            
            
            % Function 3 (columns 7-9)

            subplot(3,9, (row-1)*9 + 1 + 6);
            hold on
            plot(x, MeanGaussians{row, concentration_idx(col)}{1}, 'Color', concentration_colors(concentration_idx(col),:), 'Linewidth', 2)
            %fill([x fliplr(x)],[UpperboundGaussians{row, concentration_idx(col)}{1}, fliplr(LowerboundGaussians{row, concentration_idx(col)}{1})], ...
            %    concentration_colors(concentration_idx(col),:),'FaceAlpha', 0.2);    
            
            %stdshade(MeanGaussians{row, concentration_idx(col)}{1}, 0.2,  concentration_colors(concentration_idx(col),:), x' ) %alpha,acolor,F,smth
            subplot(3,9, (row-1)*9 + 2 + 6);
            hold on
            plot(x, MeanGaussians{row, concentration_idx(col)}{2}, 'Color', concentration_colors(concentration_idx(col),:), 'Linewidth', 2)
            %fill([x fliplr(x)],[UpperboundGaussians{row, concentration_idx(col)}{1}, fliplr(LowerboundGaussians{row, concentration_idx(col)}{1})], ...
            %    concentration_colors(concentration_idx(col),:),'FaceAlpha', 0.2);    
            subplot(3,9, (row-1)*9 + 3 + 6);
            hold on
            plot(x, MeanGaussians{row, concentration_idx(col)}{3}, 'Color', concentration_colors(concentration_idx(col),:), 'Linewidth', 2)
            %fill([x fliplr(x)],[UpperboundGaussians{row, concentration_idx(col)}{1}, fliplr(LowerboundGaussians{row, concentration_idx(col)}{1})], ...
            %    concentration_colors(concentration_idx(col),:),'FaceAlpha', 0.2);    
            
            title(sprintf('Odor %d, Conc %d', row, col));
        end
    end
    disp(fullfile('/Users/barrab01/Documents/Repos/ephysMH/figures/BoldingFigure/', ['M', num2str(imouse), '.png']))
    saveas(gcf,fullfile('/Users/barrab01/Documents/Repos/ephysMH/figures/BoldingFigure/', ['M', num2str(imouse), 'PCx_summary.png']))
end




%% USELESS

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



