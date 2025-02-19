addpath(genpath('/Users/barrab01/Documents/Repos/matnwb'));
addpath(genpath('/Users/barrab01/Documents/Repos/ephysMH'));
addpath(genpath('/Users/barrab01/Documents/Repos/chronux_2_12'));

filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_12/040136__24_11_12.nwb'
filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040136/24_11_11/040136__24_11_11.nwb'
filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/24_11_18/040118__24_11_18.nwb'
nwbFile = nwbRead(filename);
analogFS = 1000; 
dt = 1/analogFS; 
%% Summary of data structure

% nwbFile.acquisition is a Set object. A set object has properties and
% methods that you can list with properties(nwbFile.acquisition) or 
% methods(nwbFile.acquisition). 
% nwbFile.acquisition contains
% - pid_data
% - respiration_data
% which you can extract by using 'get" - nwbFile.acquisition.get('pid_data')

% nwbFile.analysis contains
% - filtered respiration data [types.core.TimeSeries]
% - respiratory_features [types.hdmf_common.DynamicTable]
% - trial_breaths [types.hdmf_common.DynamicTable]
% - trial_breaths_Pre [types.hdmf_common.DynamicTable]

% nwbFile.units contains 
% - spike_times in nwbFile.units.spike_times.data(:)
% - firing_rate in nwbFile.units.vectordata.get('firing_rate').data(:) and
% much else 

%% Small initial script 
figure
% Retrieve start and stop times
trials_start_stop_time = [nwbFile.intervals.get('odor_trials').start_time.data(:), ...
                        nwbFile.intervals.get('odor_trials').stop_time.data(:)]

Odors = nwbFile.intervals.get('odor_trials').vectordata.get('pins').data(:); 
Concentration = nwbFile.intervals.get('odor_trials').vectordata.get('concentration').data(:); 
Sniff = nwbFile.acquisition.get('respiration_data').data(:);  
SniffTime = nwbFile.acquisition.get('respiration_data').timestamps(:);  
SniffFS = 1/mean(diff(nwbFile.acquisition.get('respiration_data').timestamps(:))); 
Sniffdt = mean(diff(nwbFile.acquisition.get('respiration_data').timestamps(:)))
plot_time = 1; 
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


%% Some in-depth sniff analysis

Sniff = nwbFile.acquisition.get('respiration_data').data(:);  
SniffTime = nwbFile.acquisition.get('respiration_data').timestamps(:);  
fSniff = lowPass(Sniff, round(1/Sniffdt), 10, 3); 
baselineSniff = lowPass(Sniff, round(1/Sniffdt),0.8, 2); 

figure
plot(SniffTime, (Sniff), 'k')
hold on
plot(SniffTime(1:end-1), dSniff, 'b', 'Linewidth', 2)
plot(SniffTime,fSniff, 'r')
plot(SniffTime,baselineSniff, 'g', 'Linewidth', 2)

xlim([2510, 2520])
%% Add sniff to dataset
plot_time = 1; 
for itrial = 1 : size(trials_start_stop_time, 1)
    itrial% find sniff snip 
    [minv, idxSniff] = min(abs(SniffTime - trials_start_stop_time(itrial, 1))); 
    SniffSnip = Sniff(idxSniff : idxSniff + floor(plot_time/ Sniffdt)); 
    fSniffSnip = fSniff(idxSniff : idxSniff + floor(plot_time/ Sniffdt));
    basSniffSnip = baselineSniff(idxSniff : idxSniff + floor(plot_time/ Sniffdt));
    SniffTimeSnip = SniffTime(idxSniff : idxSniff + floor(plot_time/ Sniffdt)); 
    spiketime_trains{itrial}.Sniff = SniffSnip; 
    spiketime_trains{itrial}.fSniff = fSniffSnip; 
    spiketime_trains{itrial}.bSniff = basSniffSnip; 
    spiketime_trains{itrial}.SniffTime = SniffTimeSnip; 
end

%% Find inhalation
skipped_samples = 50; 
figure
discarded_trials = []; 
for itrial = 1 : length(spiketime_trains)
    
    % Filter Sniff 

    filtered_sniff = spiketime_trains{itrial}.fSniff;  
    hold on
    %filtered_sniff(1:50) = mean(filtered_sniff); 
    %filtered_sniff(1:100) = movmean(filtered_sniff(1:100), 10); 
    % Threshold crossing 
    %th = mean(spiketime_trains{itrial}.Sniff) - 0.5*std(spiketime_trains{itrial}.Sniff); 
    
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
    
    
    
%     
%     %
%     [B,I] = sort(abs(filtered_sniff(skipped_samples:end) - th)); % I are the indexes of all threshld crossing
%     I = I + skipped_samples; 
%     % Derivative 
%     dSniff = diff(filtered_sniff);
%     ddSniff = diff(dSniff); 
%     dddSniff = diff(ddSniff); 
%     
%     [peaks, peaks_idx] = findpeaks(-dSniff, 'MinPeakHeight', mean(-dSniff)+0.5*std(-dSniff)); 
% 
%     % Find threshold crossing that is the closest to first derivative peak 
%     [minv, minidx] = min(abs(I(1:30) - peaks_idx(1))); 
%     times_th_cross = spiketime_trains{itrial}.SniffTime(I(minidx)); 
%     mintime = spiketime_trains{itrial}.SniffTime(peaks_idx(1));
%     spiketime_trains{itrial}.inhalation_idx = I(minidx); 
%     spiketime_trains{itrial}.inhalation_time = times_th_cross; 
        % Plot
        plot(spiketime_trains{itrial}.SniffTime,filtered_sniff, 'k'  )
        plot(spiketime_trains{itrial}.SniffTime, logicSniff, 'm' )
        plot([min(spiketime_trains{itrial}.SniffTime), max(spiketime_trains{itrial}.SniffTime)], [th, th], 'r')
        plot(spiketime_trains{itrial}.inhalation_time, mean(spiketime_trains{itrial}.bSniff), 'ok', 'MarkerFaceColor','k')
        %plot([min(spiketime_trains{itrial}.SniffTime), max(spiketime_trains{itrial}.SniffTime)], [mean(spiketime_trains{itrial}.bSniff), mean(spiketime_trains{itrial}.bSniff)], 'g')
        %plot(spiketime_trains{itrial}.inhalation_time, th, 'ok', 'MarkerFaceColor','k')
        %plot([mintime, mintime], [-0.2,2], 'k--')
        %plot(spiketime_trains{itrial}.SniffTime(2:end), dSniff*10, 'g')
        %plot(spiketime_trains{itrial}.SniffTime(3:end), ddSniff*1000, 'b')
        %plot(spiketime_trains{itrial}.SniffTime(4:end), dddSniff*10000, 'r')
        %plot([spiketime_trains{itrial}.SniffTime(1), spiketime_trains{itrial}.SniffTime(end)], [mean(spiketime_trains{itrial}.Sniff), mean(spiketime_trains{itrial}.Sniff)], 'b--')
        title(itrial)
        pause(1)
        clf
    catch
        title('Trial without sniff')
        discarded_trials = [discarded_trials; itrial]; 
        plot(spiketime_trains{itrial}.SniffTime,filtered_sniff, 'k'  )
        plot(spiketime_trains{itrial}.SniffTime, logicSniff, 'm' )
        pause(1)
        clf
    end
    
    %for i = I(1:20)'
    %    plot(spiketime_trains{itrial}.SniffTime(i), th, 'o')
    %end
    
end

%% Find trials for each condition
trials_numbers= struct(); 
% vial 5 - Ethyl Tiglate 
% vial 6 - Ethyl Butyrate 
% vial 7 - Acetophenone
% vial 8 - 2-Heptanone

odornums =[5, 6, 7, 8]; 
odorlabels = {'ethyltiglate', 'ethylbutyrate', 'acetophenone', 'heptanone'}; 
for odor = unique(Odors)'
    idx = find(odornums == odor); 
    trials_numbers.(odorlabels{idx}) = struct(); 
    idx_odor = find(Odors ==odor);
    all_conc = unique(Concentration(idx_odor)); 
    trials_numbers.(odorlabels{idx}).concs = all_conc; 
    trials_numbers.(odorlabels{idx}).trial_idxs = []; 
    for iconc  = 1: length(all_conc)
        trials_numbers.(odorlabels{idx}).trial_idxs{iconc} = intersect(find(Concentration == all_conc(iconc)),idx_odor );
         
    end

    
end

%% PCA and clustering based on FR per single condition
sig = 0.020; 
mint = -0.01; 
maxt = 0.5; 
% Extract good PC and OB units
goodPCunits = strfind(spiketime_trains{itrial}.unittype, '2');
goodOBunits = strfind(spiketime_trains{itrial}.unittype, '1');
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

OBMAT = []; 
for iunit = goodOBunits
    OBline = []; 
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
            OBline= [OBline,meanFR ]; 
        end
    end
    OBMAT= [OBMAT; OBline]; 
end
OBMATraw = OBMAT; 
OBMAT = (zscore(OBMATraw')'); 

% Plot PCA and kmeans analysis for PC and OB units
clustercols = {'r', 'g', 'b', 'y'}; 
[clusterPCCells, PCIDX] = cluster_units(PCMAT, clustercols); 
[clusterOBCells, OBIDX] = cluster_units(OBMAT, clustercols); 

for i = unique(OBIDX)'
    OBClSize(i) = length(find(OBIDX == i)); 
end
[U, exOBCL] = max(OBClSize); 
for i = unique(PCIDX)'
    PCClSize(i) = length(find(PCIDX ==i)); 
end
[U, exPCCL] = max(PCClSize); 

%% Now, for the units that we have found to be part of a specific cluster (the more numerous is the one where cells
% have firing rate increasing for concentration, now re-do the histograms

%% Bulb first
sig = 0.020; 
mint = -0.2; 
maxt = 0.5; 
figOB = figure(); 
for iodor = 1: length(odorlabels)
    for iconc = 1 : length(trials_numbers.(odorlabels{iodor}).concs)
        
        allunitsPC = [];
        allunitsOB = []; 
        for ii = 1 : length(trials_numbers.(odorlabels{iodor}).trial_idxs{iconc})
            itrial = trials_numbers.(odorlabels{iodor}).trial_idxs{iconc}(ii); 
            for iunit = goodOBunits(OBIDX == exOBCL)
                idx = intersect(find(spiketime_trains{itrial}.units{iunit}> spiketime_trains{itrial}.inhalation_time+mint), ...
                    find(spiketime_trains{itrial}.units{iunit}< spiketime_trains{itrial}.inhalation_time + maxt)); 
                allunitsOB = [allunitsOB; spiketime_trains{itrial}.units{iunit}(idx)]; 
                
            end
            figure(figOB)
            subplot(length(odorlabels), length(trials_numbers.(odorlabels{iodor}).concs), (iodor-1)*length(trials_numbers.(odorlabels{iodor}).concs) + iconc)
            bins = linspace(mint, maxt, 50); 
            %[N,X] = hist(allunitsOB-spiketime_trains{itrial}.inhalation_time,bins ); 
            hold on
            histogram(allunitsOB-spiketime_trains{itrial}.inhalation_time,bins,'facealpha',.7,'edgecolor','none')
        end
        title([num2str(spiketime_trains{itrial}.odor), num2str(spiketime_trains{itrial}.conc)])
        ylim([0, 80])
        plot([0, 0], [0, 80], 'k--', 'Linewidth', 2)
    end
end
%% Piriform cortex 
figPC = figure(); 
for iodor = 1: length(odorlabels)
    for iconc = 1 : length(trials_numbers.(odorlabels{iodor}).concs)
        
        allunitsPC = [];
        for ii = 1 : length(trials_numbers.(odorlabels{iodor}).trial_idxs{iconc})
            itrial = trials_numbers.(odorlabels{iodor}).trial_idxs{iconc}(ii); 
            for iunit = goodPCunits(PCIDX == exPCCL)
                idx = intersect(find(spiketime_trains{itrial}.units{iunit}> spiketime_trains{itrial}.inhalation_time+mint), ...
                    find(spiketime_trains{itrial}.units{iunit}< spiketime_trains{itrial}.inhalation_time + maxt)); 
                allunitsPC = [allunitsPC; spiketime_trains{itrial}.units{iunit}(idx)]; 
                
            end
            figure(figPC)
            subplot(length(odorlabels), length(trials_numbers.(odorlabels{iodor}).concs), (iodor-1)*length(trials_numbers.(odorlabels{iodor}).concs) + iconc)
            bins = linspace(mint, maxt, 50); 
            %[N,X] = hist(allunitsOB-spiketime_trains{itrial}.inhalation_time,bins ); 
            hold on
            histogram(allunitsPC-spiketime_trains{itrial}.inhalation_time,bins,'facealpha',.7,'edgecolor','none')
        end
        title([num2str(spiketime_trains{itrial}.odor), num2str(spiketime_trains{itrial}.conc)])
        ylim([0, 80])
        plot([0, 0], [0, 80], 'k--', 'Linewidth', 2)
    end
end


%% Looking at Sniff 
for iodor = 1: length(odorlabels)
    figmean = figure(); 
    figPC = figure('units','normalized','outerposition',[0 0 1 0.5])
    hold on
    for iconc = 1 : length(trials_numbers.(odorlabels{iodor}).concs)
        subplot(1, length(trials_numbers.(odorlabels{iodor}).trial_idxs{iconc}), iconc)
        hold on
        for ii = 1 : length(trials_numbers.(odorlabels{iodor}).trial_idxs{iconc})
            
            itrial = trials_numbers.(odorlabels{iodor}).trial_idxs{iconc}(ii);
            times = spiketime_trains{itrial}.SniffTime(spiketime_trains{itrial}.inhalation_idx:spiketime_trains{itrial}.inhalation_idx+300) - spiketime_trains{itrial}.SniffTime(spiketime_trains{itrial}.inhalation_idx)
            sniff =  spiketime_trains{itrial}.Sniff(spiketime_trains{itrial}.inhalation_idx:spiketime_trains{itrial}.inhalation_idx+300)
            plot(times, sniff, 'k')
            
        
        end
        
    end
end

%%
k = 3; 
c = hot(12);


for iodor = 1: length(odorlabels)
    figmean = figure(); 
    figPC = figure('units','normalized','outerposition',[0 0 1 0.5])
    hold on
    for iconc = 1 : length(trials_numbers.(odorlabels{iodor}).concs)
        
        MU = []; 
        SIGMA = []; 
        for ii = 1 : length(trials_numbers.(odorlabels{iodor}).trial_idxs{iconc})
            allunitsPC = [];
            itrial = trials_numbers.(odorlabels{iodor}).trial_idxs{iconc}(ii); 
            for iunit = goodPCunits(PCIDX == exPCCL)
                idx = intersect(find(spiketime_trains{itrial}.units{iunit}> spiketime_trains{itrial}.inhalation_time+mint), ...
                    find(spiketime_trains{itrial}.units{iunit}< spiketime_trains{itrial}.inhalation_time + maxt)); 
                allunitsPC = [allunitsPC; spiketime_trains{itrial}.units{iunit}(idx)]; 
                
            end
        
            X = allunitsPC-(spiketime_trains{itrial}.inhalation_time); 
            [mu, sigma] = GMM_clustering(X,k); 
            disp([mu, sigma])
            [mu, idx]=sort(mu); 
            sigma = sigma(idx)
            disp([mu, sigma])
            MU = [MU; mu]; 
            SIGMA = [SIGMA; sigma]; 
            % Plotting 
            figure(figPC)
            subplot(1, length(trials_numbers.(odorlabels{iodor}).trial_idxs{iconc}), ii)
            plot(spiketime_trains{itrial}.SniffTime - spiketime_trains{itrial}.inhalation_time, 100*spiketime_trains{itrial}.Sniff)
            hold on
            %hist(X, 100)
            %hold on
            x = [0:0.001:1];
            % Plot over the existing figure, using black lines for the estimated pdfs.
            for i = 1 : k
        
                y1 = gaussian1D(x, mu(i), sigma(i));
                plot(x, y1, '-','Color', c(iconc,:) , 'Linewidth', 0.5);
            end
            ylim([0, 50])
        end
        
        figure(figmean)
        hold on
        for i = 1 : 3
            y1 = gaussian1D(x, mean(MU(:,i)), mean(SIGMA(:, i)));
            plot(x, y1, '-','Color', c(iconc,:) , 'Linewidth', 2);
        end
        title([num2str(odorlabels{iodor}), '_' , num2str(trials_numbers.(odorlabels{iodor}).concs(iconc))])
                
        
    end
    pause
    
end
%% Working out histogram fitting
k = 3;
figure
for itrial = 1:290
    
  
    allunitsPC = []; 
    mint = 0; 
    maxt = 0.6; 
    for iunit = goodPCunits(PCIDX == exPCCL)
        idx = intersect(find(spiketime_trains{itrial}.units{iunit}> spiketime_trains{itrial}.inhalation_time+mint), ...
            find(spiketime_trains{itrial}.units{iunit}< spiketime_trains{itrial}.inhalation_time + maxt)); 
        allunitsPC = [allunitsPC; spiketime_trains{itrial}.units{iunit}(idx)]; 
        
    end
    %histfit(allunitsPC-(spiketime_trains{itrial}.inhalation_time ),50, 'weibull')
   
    X = allunitsPC-(spiketime_trains{itrial}.inhalation_time); 
    % Set 'k' to the number of clusters to find.
    
    [mu, sigma] = GMM_clustering(X,k); 
    %=====================================================
    hist(X, 100)
    hold on
    x = [0:0.001:1];
    
   
    % Plot over the existing figure, using black lines for the estimated pdfs.
    for i = 1 : k

        y1 = gaussian1D(x, mu(i), sigma(i));
        plot(x, y1, 'r-', 'Linewidth', 2);
    end

    pause
    clf
end


%% trying gaussian mixture models 
% Set 'm' to the number of data points.
X = allunitsPC-(spiketime_trains{itrial}.inhalation_time); 
% Set 'k' to the number of clusters to find.

[mu, sigma] = GMM_clustering(X,k); 
%=====================================================
hist(X, 100)
hold on
x = [0:0.001:0.6];
y1 = gaussian1D(x, mu(1), sigma(1));
y2 = gaussian1D(x, mu(2), sigma(2));

% Plot over the existing figure, using black lines for the estimated pdfs.
plot(x, y1, 'k-');
plot(x, y2, 'k-');



%% Plot single-cell, single-condition raster - this is the bolding plot

iunit = 1206; 
iodor = 1; 
iconc = 3;
sig = 0.020; 
mint = -0.2; 
maxt = 0.5; 
figure
goodPCunits = strfind(spiketime_trains{itrial}.unittype, '2');
goodOBunits = strfind(spiketime_trains{itrial}.unittype, '1');
for iunit = goodPCunits
    for iodor = 1: length(odorlabels)
        for iconc = 1 : length(trials_numbers.(odorlabels{iodor}).concs)
            subplot(length(odorlabels), length(trials_numbers.(odorlabels{iodor}).concs), (iodor-1)*length(trials_numbers.(odorlabels{iodor}).concs) + iconc)
            try
                Rast = {}; 
                for ii = 1 : length(trials_numbers.(odorlabels{iodor}).trial_idxs{iconc})
                    itrial = trials_numbers.(odorlabels{iodor}).trial_idxs{iconc}(ii); 
                    idx = intersect(find(spiketime_trains{itrial}.units{iunit}> spiketime_trains{itrial}.inhalation_time-0.2), ...
                                        find(spiketime_trains{itrial}.units{iunit}< spiketime_trains{itrial}.inhalation_time + 0.5)); 
                    spikes = spiketime_trains{itrial}.units{iunit}(idx)-spiketime_trains{itrial}.inhalation_time; 
                    Rast{ii} = spikes; 
                    for spike = spikes'
                        plot([spike, spike ], [-ii, -ii-0.8], 'k')
                        hold on
                    end
                
                end
                [iunit, iodor, iconc]
                tvect = linspace(mint, maxt, (maxt-mint)/sig); 
                [FR,t, err] = compute_FR_from_raster(Rast , sig, mint, maxt, tvect ); 
                plot(t, FR, 'b')
                plot([0, 0], [0, 10], 'r--')
                ylim([-20, 60])
            catch
                disp(['problem with unit'])
                title(['PROB'])
            end
        end
    end
    pause
    clf

end


%% Plot al trials mixed together (not sure this is correct)

figPC = figure(); 
figOB = figure(); 
for iodor = 1: length(odorlabels)
    for iconc = 1 : length(trials_numbers.(odorlabels{iodor}).concs)
        itrial = trials_numbers.(odorlabels{iodor}).trial_idxs{iconc}(1); 
        allunitsPC = [];
        allunitsOB = []; 
        for ii = 1 : length(trials_numbers.(odorlabels{iodor}).trial_idxs{iconc})
            itrial = trials_numbers.(odorlabels{iodor}).trial_idxs{iconc}(ii); 
            for iunit = 1 : length(spiketime_trains{itrial}.units)
                if strcmp(spiketime_trains{itrial}.unittype(iunit), '2')
                    idx = intersect(find(spiketime_trains{itrial}.units{iunit}> spiketime_trains{itrial}.inhalation_time), ...
                        find(spiketime_trains{itrial}.units{iunit}< spiketime_trains{itrial}.inhalation_time + 0.2)); 
                    allunitsPC = [allunitsPC; spiketime_trains{itrial}.units{iunit}(idx)-spiketime_trains{itrial}.inhalation_time]; 
                end
                if strcmp(spiketime_trains{itrial}.unittype(iunit), '1')
                    idx = intersect(find(spiketime_trains{itrial}.units{iunit}> spiketime_trains{itrial}.inhalation_time), ...
                        find(spiketime_trains{itrial}.units{iunit}< spiketime_trains{itrial}.inhalation_time + 0.2)); 
                    allunitsOB = [allunitsOB; spiketime_trains{itrial}.units{iunit}(idx)-spiketime_trains{itrial}.inhalation_time]; 
                end
            end
        end
        figure(figOB)
        bins = linspace(0, 0.2, 50); 
        subplot(4,8, (iodor-1)*length(all_conc) + iconc)
        hist(allunitsOB,bins )
        title([odorlabels{iodor}, '  ', num2str(all_conc(iconc))])
        ylim([0, 120])
        figure(figPC)
        subplot(4,8, (iodor-1)*length(all_conc) + iconc)        
        hist(allunitsPC, bins)
        title([odorlabels{iodor}, '  ', num2str(all_conc(iconc))])
        ylim([0, 120])
    end
end


%% Look at repeatibility of conditions

figOB = figure(); 
figPC = figure(); 
figSN = figure(); 
iodor = 2; 
for iconc = 1 : length(trials_numbers.(odorlabels{iodor}).concs)
    itrial = trials_numbers.(odorlabels{iodor}).trial_idxs{iconc}(1); 
    
    for ii = 1 : length(trials_numbers.(odorlabels{iodor}).trial_idxs{iconc})
        allunitsPC = [];
        allunitsOB = []; 
        itrial = trials_numbers.(odorlabels{iodor}).trial_idxs{iconc}(ii); 
        for iunit = 1 : length(spiketime_trains{itrial}.units)
            if strcmp(spiketime_trains{itrial}.unittype(iunit), '2')
                idx = intersect(find(spiketime_trains{itrial}.units{iunit}> spiketime_trains{itrial}.inhalation_time), ...
                    find(spiketime_trains{itrial}.units{iunit}< spiketime_trains{itrial}.inhalation_time + 0.2)); 
                allunitsPC = [allunitsPC; spiketime_trains{itrial}.units{iunit}(idx)-spiketime_trains{itrial}.inhalation_time]; 
            end
            if strcmp(spiketime_trains{itrial}.unittype(iunit), '1')
                idx = intersect(find(spiketime_trains{itrial}.units{iunit}> spiketime_trains{itrial}.inhalation_time), ...
                    find(spiketime_trains{itrial}.units{iunit}< spiketime_trains{itrial}.inhalation_time + 0.2)); 
                allunitsOB = [allunitsOB; spiketime_trains{itrial}.units{iunit}(idx)-spiketime_trains{itrial}.inhalation_time]; 
            end
        end
    
        figure(figOB)
        bins = linspace(0, 0.2, 50); 
        subplot(10,8, (ii-1)*length(all_conc) + iconc)
        hist(allunitsOB,bins )
        title([odorlabels{iodor}, '  ', num2str(all_conc(iconc))])
        ylim([0, 30])
        figure(figPC)
        subplot(10,8, (ii-1)*length(all_conc) + iconc)        
        hist(allunitsPC, bins)
        hold on
        title([odorlabels{iodor}, '  ', num2str(all_conc(iconc))])
        [V, idx_inhalation] = min(abs(spiketime_trains{itrial}.SniffTime - spiketime_trains{itrial}.inhalation_time)); 
        plot(spiketime_trains{itrial}.SniffTime(idx_inhalation : idx_inhalation+200)-spiketime_trains{itrial}.SniffTime(idx_inhalation), 20+(10*spiketime_trains{itrial}.Sniff(idx_inhalation : idx_inhalation+200)))
        %ylim([0, 30])
        %figure(figSN)
        %subplot(10,8, (ii-1)*length(all_conc) + iconc)      
%         [V, idx_inhalation] = min(abs(spiketime_trains{itrial}.SniffTime - spiketime_trains{itrial}.inhalation_time)); 
%         plot(spiketime_trains{itrial}.SniffTime(idx_inhalation : idx_inhalation+200), spiketime_trains{itrial}.Sniff(idx_inhalation : idx_inhalation+200))
%         ylim([0.5, 2])
    end
end

%% Plot PSTH for different conditions, all units mixed

% YOU HAVE TO FILTER PC from OB units
figPC = figure(); 
figOB = figure(); 
for itrial = 1 : length(spiketime_trains)
    allunitsPC = [];
    allunitsOB = []; 
    for iunit = 1 : length(spiketime_trains{itrial}.units)
        if strcmp(spiketime_trains{itrial}.unittype(iunit), '2')
            idx = intersect(find(spiketime_trains{itrial}.units{iunit}> spiketime_trains{itrial}.inhalation_time), ...
                find(spiketime_trains{itrial}.units{iunit}< spiketime_trains{itrial}.inhalation_time + 0.2)); 
            allunitsPC = [allunitsPC; spiketime_trains{itrial}.units{iunit}(idx)]; 
        end
        if strcmp(spiketime_trains{itrial}.unittype(iunit), '1')
            idx = intersect(find(spiketime_trains{itrial}.units{iunit}> spiketime_trains{itrial}.inhalation_time), ...
                find(spiketime_trains{itrial}.units{iunit}< spiketime_trains{itrial}.inhalation_time + 0.2)); 
            allunitsOB = [allunitsOB; spiketime_trains{itrial}.units{iunit}(idx)]; 
        end
    end
    figure(figOB)
    bins = linspace(0, 0.2, 50); 
    subplot(10,10, itrial)
    hist(allunitsOB-spiketime_trains{itrial}.inhalation_time,bins )
    title([num2str(spiketime_trains{itrial}.odor), num2str(spiketime_trains{itrial}.conc)])
    figure(figPC)
    subplot(10,10, itrial)
    hist(allunitsPC-spiketime_trains{itrial}.inhalation_time, bins)
    title([num2str(spiketime_trains{itrial}.odor), num2str(spiketime_trains{itrial}.conc)])
end

%% if histograms are similar, are units doing the same thing from trial to trial? 

%% Do histograms suggest synchrony?

%% ISI patterns? 

%% Cross correlation across different units?? 

%% Decoder on histograms? Spike trains? Spike trains with unit id assigned? Do they overperform spike trains with no unit assignment?

%%
format long
previous_unit_idx = 1; 
for iunit = 1 : 5%max(nwbFile.units.vectordata.get('cluster_id').data(:))-1 
    thisunit_idx =  previous_unit_idx; 
    nextunit_idx =  nwbFile.units.spike_times_index.data(iunit); 
    allspikes_unit = nwbFile.units.spike_times.data(thisunit_idx: nextunit_idx); 
    [allspikes_unit(1), allspikes_unit(end)]
    previous_unit_idx = nextunit_idx+1; 
end
%% 

nwbFile.units.spike_times
nwbFile.units.vectordata

size(nwbFile.units.vectordata.get('cluster_id').data(:))

plot(nwbFile.units.vectordata.get('firing_rate').data(:))
disp(nwbFile.units.spike_times.data(:))

%%
size(nwbFile.units.vectordata.get('cell_type').data(:))
nwbFile.units.vectordata


%%

nwbFile.analysis.get('trial_breaths').vectordata.resp_features_id(:)




%%
Sniff = nwbFile.acquisition.get('respiration_data').data;  
PID = nwbFile.acquisition.get('pid_data').data;  




