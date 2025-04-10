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



%% Read the dataset, extract variables and perform quality check
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



%% Processs all mice together

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


%% Load already processed data
savefilenameallmice = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/all_mice_dataset.mat'; 

load(savefilenameallmice)

sig = 0.010; 
mint = 0; 
maxt = 0.5;
tvect = linspace(mint, maxt, (maxt-mint)/sig); 

%% Computing units that are significant
Nconditions = 32; 
clear unit_significance
figure
for iunit = 1 : size(peak_latency, 1)
    baseline_rates = cell(1, length(unique(all_mice_condition_table(:,1)))); 
    stimulus_rates = cell(1, length(unique(all_mice_condition_table(:,1)))); 
    for icondition = 1: Nconditions
        iodor = all_mice_condition_table(icondition, 1); 
        iconc = all_mice_condition_table(icondition, 2); 
        
        for itrial = 1 : length(all_mice_spikes_table{iunit, icondition})
            
            num_spikes = length(intersect(find(all_mice_spikes_table{iunit, icondition}{itrial}>0), find(all_mice_spikes_table{iunit, icondition}{itrial}<0.5))); 
            if iconc == 1
                baseline_rates{iodor} = [baseline_rates{iodor}, num_spikes]; 
            else 
                stimulus_rates{iodor} = [stimulus_rates{iodor}, num_spikes]; 
            end

        end
    end

    for iodor = 1 : length(unique(all_mice_condition_table(:, 1)))
        labels = [repmat([0], 1, length(baseline_rates{iodor})) ,  repmat([1], 1,length(stimulus_rates{iodor}))];  
        scores = [baseline_rates{iodor} , stimulus_rates{iodor}];   
        resp_sign_matrix(iunit, iodor) = sign(mean(stimulus_rates{iodor})-mean(baseline_rates{iodor})); 
        [X, Y, AUC, p] = ranksumROC(labels,scores);
        hold on
        plot(X,Y)
        unit_significance(iunit, iodor) = p; 
    end
end
figure
imagesc(unit_significance, [0, 0.05])

%% Classify units depending on responses

figure
for j = 1 : size(unit_significance, 2)
    subplot(1, size(unit_significance, 2), j)
    hist(unit_significance(:, j), 100)
    hold on
    plot([0.05, 0.05], [0, 100], 'r--')
    xlim([-0.1, 1])

    idx_unit_modulated(:,j)= (unit_significance(:, j)< 0.05); 

end

% divide in 1) units that are never modulated 2) units modulated by all 3)
% units modulated by only one odor 3) units modulated by more than one but
% less than all


units_modulated_by_one = any(idx_unit_modulated, 2); 
units_modulated_by_all = all(idx_unit_modulated, 2); 
nodors_modulating = sum(idx_unit_modulated, 2);  

units_never_modulated = length(find(units_modulated_by_one ==0)); 
units_always_modulated = length(find(units_modulated_by_all ==1)); 
units_modulated_by2odorsmin = length(intersect(find(nodors_modulating>1), find(nodors_modulating<4))); 

[units_never_modulated, units_always_modulated, sum(units_modulated_by_one)-units_always_modulated-units_modulated_by2odorsmin, units_modulated_by2odorsmin, size(idx_unit_modulated, 1)]
figure
ax(1) = subplot(1, 2, 1)
pie([units_never_modulated, units_always_modulated, sum(units_modulated_by_one)-units_always_modulated-units_modulated_by2odorsmin, units_modulated_by2odorsmin,])
labels = {'NO RESP', 'ALL', 'ONE ODOR', 'TWO/THREE ODORS'}
legend(labels)
title('Do units repond?')
colormap(ax(1) ,'hot')
% Now, of all units that respond to at least one odor, how many positive
% and negative
unit_odor_sign = sum(resp_sign_matrix(find(units_modulated_by_one), :), 2)

NEG = length(find(unit_odor_sign ==-4)); 
MNEG = length(intersect(find(unit_odor_sign >-4), find(unit_odor_sign <0))); 
NEU = length(find(unit_odor_sign ==0)); 
MPOS = length(intersect(find(unit_odor_sign >0), find(unit_odor_sign <4))); 
POS = length(find(unit_odor_sign ==4)); 
ax(2) = subplot(1, 2, 2)

pie([NEG, MNEG , NEU , MPOS , POS])% , length(unit_odor_sign)]
labels = {'NEG', 'MNEG', 'NEU', 'MPOS', 'POS'}; 
legend(labels, 'Location', 'South')
mymap = blue_white_red(5); 
colormap(ax(2) ,mymap)
title('Response sign')


% Create list of positive responding units 
sum(resp_sign_matrix, 2)
posunits = intersect(find(sum(resp_sign_matrix, 2) >0), find(sum(resp_sign_matrix, 2) <4)); 
posunits = intersect(posunits, find(units_modulated_by_one)); 
posunits


%% Plot cell specific PSTH

figure
icondition = 3; 
iodor = 1; 
[sortedpeaks, sorted_idx_this_condition] = sort(peak_latency(:, icondition)); 

Nunits = 30; 
Nconditions = 8; 
unit_multiplier = 100; 
iunit = 0; 
iiunit = 0; 
iiiunit = 0
%for iunit = 1:size(peak_latency, 1)%sorted_idx_this_condition((unit_multiplier-1)*Nunits + 1:unit_multiplier* Nunits)'% (unit_multiplier-1)*Nunits + 1 :unit_multiplier*Nunits %
while(iunit<Nunits)
    iiiunit = iiiunit + 1; 
    % MESS WITH IUNIT AND IIUNIT THINK ABOUT IT 
    iiunit = posunits(iiiunit); 
    baseline_rates = cell(length(unique(all_mice_condition_table(:,1)))); 
    stimulus_rates = cell(length(unique(all_mice_condition_table(:,1)))); 
    iicondition = 0; 
    
    for icondition = 1: Nconditions
        iodor = all_mice_condition_table(icondition, 1); 
        iconc = all_mice_condition_table(icondition, 2); 
        

        % Include only if unit responds significantly different than
        % baseline 
        isunitin= intersect(find(unit_significance(:, iodor)< 0.05), iiunit); 
        
        if ~isempty(isunitin)
            if iconc == 1
                f = figure
                set(f, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Fullscreen

                iunit = iunit +1; 
            end
            
            iicondition = iicondition +1; 
            %subplot(10, 8, (iunit-1)*8 + icondition)
            
            subplot( 1, Nconditions,(iicondition))
            hold on
            plot(tvect, allmice_FRmat{icondition}(iiunit, :), 'b')
            if icondition == 5
                plot([peak_latency(iiunit, icondition), peak_latency(iiunit, icondition)], [0, 50], 'r--')
            end
            for itrial = 1 : length(all_mice_spikes_table{iiunit, icondition})% number of trials
                spikesinwindow = all_mice_spikes_table{iiunit, icondition}{itrial}(intersect(find(all_mice_spikes_table{iiunit, icondition}{itrial}<maxt), find(all_mice_spikes_table{iiunit, icondition}{itrial}>mint))); 
                plot(spikesinwindow, -3*itrial*ones(size(spikesinwindow/sig)), '|k', 'MarkerSize',2)
            end
            ylim([-50, 50])
    %         if icondition == 5
    %              xlim([0, 0.05])
    %              ylim([0, max(allmice_FRmat{icondition}(iunit, :)+2)])
    %         end
            title(['U ', num2str(iiunit), ',_ S', num2str(icondition)])
        
        end
    end
end



%% Plot PSTHs for all units , all mice 
clear sorted_responses
Nodors = 4
Nconcs =8; 
odorcols = ['r', 'b', 'y', 'g']; 
figureFR = figure()
figureH= figure()
%figureC= figure()
black2bluemap = blackToBlue(10); 
black2redmap = blackToRed(10); 
yellow2redmap = yellowToRed(10); 
odor_onset = abs(mint/sig); 
ms100 = odor_onset + abs(0.1/sig); 
ms200 = odor_onset + abs(0.2/sig); 
    
%Sort responses across all mice
% iconc = 5; 
% for iodor = 1 : Nodors % for each odor
%     idxcondition = intersect(find(all_mice_condition_table(:,1) == iodor), find(all_mice_condition_table(:,2) == iconc)); 
%     % sort for each condition separately
%     icondition = mean(all_mice_condition_table(idxcondition, 2)); 
%     [sortedpeaks, sorted_idx{iodor}] = sort(peak_latency(:, icondition)); 
% end
figtraces = figure()
clear sorted_responses
for icondition = 1 : size(peak_latency, 2)
    
  
    iodor = all_mice_condition_table(icondition, 1); 
    iconc = all_mice_condition_table(icondition, 2); 
    %This two lines above I can do just because that all mice are ordered the same and have the same conditions
    
    % sort for each condition separately
    [sortedpeaks, sorted_idx{iodor}] = sort(peak_latency(:, icondition));  % each
    ii = 0; 
    newpeaks = []; 
    for i = sorted_idx{iodor}'
        ii = ii + 1; 
        if ~isnan(peak_latency(i, icondition)) %Check that there is peak
            sorted_responses(ii, :) = allmice_FRmat{icondition}(i, :); 
            newpeaks(ii) = peak_latency(i, icondition);
        end
    end
    
    % Firing rates curves
    %figure(figureC)
    %subplot(Nodors+2, Nconcs, (iodor-1)*Nconcs + iconc)
    %plot(sorted_responses')

    % Firing rates
    figure(figureFR)
    subplot(Nodors+2, Nconcs, (iodor-1)*Nconcs + iconc)
    imagesc((sorted_responses), [0,30])
    colormap(flipud(pink))
    
    hold on
    idx_inh = floor(size(sorted_responses, 2)*abs(-mint/(-mint + maxt))); 
    plot([idx_inh, idx_inh], [0, size(sorted_responses, 1)], '--w', 'Linewidth', 2) 
    plot(newpeaks/sig, linspace(1, length(newpeaks), length(newpeaks)), '.', 'Color', 'black', "Markersize", 2)
    plot([odor_onset, odor_onset], [0 size(sorted_responses, 1)], 'k--')
    plot([ms100, ms100], [0 size(sorted_responses, 1)], 'k--')
    plot([ms200, ms200], [0 size(sorted_responses, 1)], 'k--')
    
    % overlap of odorants
    subplot(Nodors+2, Nconcs, (Nodors)*Nconcs + iconc)
    hold on
    plot(newpeaks/sig, fliplr(linspace(1, length(newpeaks), length(newpeaks))), 'o', 'Color', odorcols(iodor), "Markersize", 1)
    % Histograms
    figure(figureH)
    subplot(Nodors, Nconcs, (iodor-1)*Nconcs + iconc)

    hist(newpeaks, 40)
    title([num2str(iodor) , num2str(iconc) , ])
    ylim([0, 1200])

    figure(figtraces)
    subplot(1, 4, iodor)
    hold on
    if iconc ==1
         plot(mean((sorted_responses)), 'Color', 'k')
    else
        plot(mean((sorted_responses)), 'Color', yellow2redmap(iconc,:))
    end

end



%% Plot average response for units that are positive and negative separately
clear sorted_responses
Nodors = 4
Nconcs =8; 
odorcols = ['r', 'b', 'y', 'g']; 
figureFR = figure()
figureH= figure()
%figureC= figure()

black2redmap = blackToRed(10); 
yellow2redmap = yellowToRed(10); 
odor_onset = abs(mint/sig); 
ms100 = odor_onset + abs(0.1/sig); 
ms200 = odor_onset + abs(0.2/sig); 
set(gcf, 'Renderer', 'opengl');
figtraces = figure()
clear sorted_responses
for icondition = 1 : size(peak_latency, 2)
    
  
    iodor = all_mice_condition_table(icondition, 1); 
    iconc = all_mice_condition_table(icondition, 2); 
    %This two lines above I can do just because that all mice are ordered the same and have the same conditions
    
    % sort for each condition separately
    [sortedpeaks, sorted_idx{iodor}] = sort(peak_latency(:, icondition));  % each
    ii = 0; 
    newpeaks = []; 
    for i = sorted_idx{iodor}'
        ii = ii + 1; 
        if ~isnan(peak_latency(i, icondition))  && resp_sign_matrix(i, iodor) ==-1 %Check that there is peak and the cell is positive-responding
            sorted_responses(ii, :) = allmice_FRmat{icondition}(i, :); 
            newpeaks(ii) = peak_latency(i, icondition);
        end
    end
    
    % Firing rates curves
    %figure(figureC)
    %subplot(Nodors+2, Nconcs, (iodor-1)*Nconcs + iconc)
    %plot(sorted_responses')

    % Firing rates
    figure(figureFR)
    subplot(Nodors+2, Nconcs, (iodor-1)*Nconcs + iconc)
    imagesc((sorted_responses), [0,30])
    colormap(flipud(pink))
    
    hold on
    idx_inh = floor(size(sorted_responses, 2)*abs(-mint/(-mint + maxt))); 
    plot([idx_inh, idx_inh], [0, size(sorted_responses, 1)], '--w', 'Linewidth', 2) 
    scatter(newpeaks/sig, linspace(1, length(newpeaks), length(newpeaks)), 1, 'k', 'filled', 'Marker', '.', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.8)
    plot([odor_onset, odor_onset], [0 size(sorted_responses, 1)], 'k--')
    plot([ms100, ms100], [0 size(sorted_responses, 1)], 'k--')
    plot([ms200, ms200], [0 size(sorted_responses, 1)], 'k--')
    

    figure(figtraces)
    subplot(1, 4, iodor)
    hold on
    if iconc ==1
         plot(mean((sorted_responses)), 'Color', 'k')
    else
        plot(mean((sorted_responses)), 'Color', yellow2redmap(iconc,:))
    end

end

%% Gaussian Mixture analysis on latency of firing rates peaks histogram


Nodors = 4
Nconcs =8; 
odorcols = ['r', 'b', 'y', 'g']; 

figureH= figure()

clear sorted_responses
nGauss = 2; 
figoverlap = figure; 
for icondition = 1 : size(peak_latency, 2)
    
  
    iodor = all_mice_condition_table(icondition, 1); 
    iconc = all_mice_condition_table(icondition, 2); 
    
    % Filter units that are significantly different from baseline
    idx_responses  = intersect(find(~isnan(peak_latency(:, icondition))) , find( resp_sign_matrix(:, iodor) ==1)); 
    idx_responses  =find(~isnan(peak_latency(:, icondition))); 
    
    X = peak_latency(idx_responses, icondition); 
    x = [-0.1:0.001:0.5];
    xRange = x; 
    % Histograms
    figure(figureH)
    subplot(Nodors, Nconcs, (iodor-1)*Nconcs + iconc)
    hold on
    
    histogram(X, 'Normalization', 'pdf', 'BinWidth', 0.01)
    hold on
    %ksdensity(X)
    [f,xi] = ksdensity(X);
    plot(xi,f, 'Linewidth', 2, 'Color', 'b');
    options = statset('MaxIter',10000);
    gmmModel = fitgmdist(X, nGauss, 'Start', 'plus', 'CovarianceType', 'full', 'Options', options);
    pdfValues = pdf(gmmModel, xRange');
    %plot(x, pdfValues, 'r-', 'Linewidth', 2); %plotting the sum of all
    %components
    mixingProportions = gmmModel.ComponentProportion;  % Mixing proportions
    means = gmmModel.mu;  % Means of the components
    covariances = gmmModel.Sigma;
    for i = 1:length(mixingProportions)
        % Create the individual Gaussian curve using the mean and covariance
        y = mixingProportions(i) * normpdf(x, means(i), sqrt(covariances(i)));  % PDF of each Gaussian
        plot(x, y, 'r', 'LineWidth', 2);  % Plot the individual Gaussian component
    end
    [v, idx_ordered_gauss] = sort(means); 
    
%     [mu, sigma] = GMM_clustering(X,nGauss); 
%     for i = 1 : nGauss
%         y1 = gaussian1D(x, mu(i), sigma(i));
%         plot(x, y1, 'g-', 'Linewidth', 2);
%     end

    figure(figoverlap)
    index_gaussian=0; 
    for i =idx_ordered_gauss'
        index_gaussian = index_gaussian+1; 
        y1 = gaussian1D(x, means(i), sqrt(covariances(i)));
    
        subplot(3, 4, iodor)
        hold on
        if index_gaussian ==1
            plot(x, y1, '-', 'Linewidth', 2, 'Color', black2redmap(iconc,:));
        else
            plot(x, y1, '-', 'Linewidth', 2, 'Color', black2bluemap(iconc,:));
        end
        
        subplot(3, 4, 8 + iodor)
        %ecdf();
        [F,xcdf,FLO,FUP] = ecdf(X); 
        hold on
        plot(xcdf,F, 'Color', black2redmap(iconc,:))
        [maxv, maxidx ] = max(y1);
        avglatency(iodor,iconc, index_gaussian) = maxidx*mean(diff(x)); 
        
    end
    title([num2str(iodor) , num2str(iconc) , ])
    %ylim([0, 1])
end

figure(figoverlap)
markers = []
for iodor = 1: 4
    subplot(3, 4, 4 + iodor)
    hold on
    for ig = 1 : 2
        plot([1:8], avglatency(iodor,:, ig), 'ko-')
    
    end
end



%% Gaussian Mixture analysis on latency of peak of spiking histograms 


Nodors = 4
Nconcs =8; 
odorcols = ['r', 'b', 'y', 'g']; 

figureH= figure()

clear sorted_responses
nGauss = 3; 
figoverlap = figure; 
for icondition = 1 : size(peak_latency, 2)
    
  
    iodor = all_mice_condition_table(icondition, 1); 
    iconc = all_mice_condition_table(icondition, 2); 
    
    % Filter units that are significantly different from baseline
    idx_responses  = intersect(find(~isnan(peak_latency(:, icondition))) , find( resp_sign_matrix(:, iodor) ==1)); 
    idx_responses  =find(~isnan(peak_latency(:, icondition))); 
    
    % extract all spikes
    X = []; 
    for iunit =  idx_responses'
        for itrial = 1 : length(all_mice_spikes_table{iunit, icondition})
            idx = intersect(find(all_mice_spikes_table{iunit, icondition}{itrial}>0), find(all_mice_spikes_table{iunit, icondition}{itrial}<0.5)); 
            X = [X; all_mice_spikes_table{iunit, icondition}{itrial}(idx)]; 
        end
    end
    
    x = [-0.1:0.001:0.5];
    xRange = x; 
    % Histograms
    figure(figureH)
    subplot(Nodors, Nconcs, (iodor-1)*Nconcs + iconc)
    hold on
    
    histogram(X, 'Normalization', 'pdf', 'BinWidth', 0.01)
    hold on
    %ksdensity(X)
    [f,xi] = ksdensity(X);
    plot(xi,f, 'Linewidth', 2, 'Color', 'b');
    options = statset('MaxIter',10000);
    gmmModel = fitgmdist(X, nGauss, 'Start', 'plus', 'CovarianceType', 'full', 'Options', options);
    pdfValues = pdf(gmmModel, xRange');
    %plot(x, pdfValues, 'r-', 'Linewidth', 2); %plotting the sum of all
    %components
    mixingProportions = gmmModel.ComponentProportion;  % Mixing proportions
    means = gmmModel.mu;  % Means of the components
    covariances = gmmModel.Sigma;
    for i = 1:length(mixingProportions)
        % Create the individual Gaussian curve using the mean and covariance
        y = mixingProportions(i) * normpdf(x, means(i), sqrt(covariances(i)));  % PDF of each Gaussian
        plot(x, y, 'r', 'LineWidth', 2);  % Plot the individual Gaussian component
    end
    [v, idx_ordered_gauss] = sort(means); 
    
%     [mu, sigma] = GMM_clustering(X,nGauss); 
%     for i = 1 : nGauss
%         y1 = gaussian1D(x, mu(i), sigma(i));
%         plot(x, y1, 'g-', 'Linewidth', 2);
%     end

    figure(figoverlap)
    index_gaussian=0; 
    for i =idx_ordered_gauss'
        index_gaussian = index_gaussian+1; 
        y1 = gaussian1D(x, means(i), sqrt(covariances(i)));
        
        subplot(3, 4, iodor)
        hold on
        if index_gaussian ==1
            plot(x, y1, '-', 'Linewidth', 2, 'Color', black2redmap(iconc,:));
        else
            plot(x, y1, '-', 'Linewidth', 2, 'Color', black2bluemap(iconc,:));
        end
        
        subplot(3, 4, 8 + iodor)
        %ecdf();
        [F,xcdf,FLO,FUP] = ecdf(X); 
        hold on
        plot(xcdf,F, 'Color', black2redmap(iconc,:))
        [maxv, maxidx ] = max(y1);
        avglatency(iodor,iconc, index_gaussian) = maxidx*mean(diff(x)); 
        
    end

    
    title([num2str(iodor) , num2str(iconc) , ])
    %ylim([0, 1])
end

figure(figoverlap)
markers = []
for iodor = 1: 4
    subplot(3, 4, 4 + iodor)
    hold on
    for ig = 1 : 2
        plot([1:8], avglatency(iodor,:, ig), 'ko-')
    
    end
end

%% Intensity decoder from gaussians built on spike histograms
%  Different trials are separated into different data points
%
% Read the spikes features
x = [-0.1:0.001:0.5];
xRange = x; 
nGauss = 3; 
gfeats = []; 
conditions = [2:8, 10:16, 18:24, 26:32]; % skipping no stim
concentration_label = []; 
for icondition = conditions
    % Find parameters
    iodor = all_mice_condition_table(icondition, 1); 
    iconc = all_mice_condition_table(icondition, 2); 
    concentration = all_mice_condition_table(icondition, 3); 
    
    % Filter units that have a peak
    idx_responses  =find(~isnan(peak_latency(:, icondition))); 
    
    % Extract all spikes
    X = repmat({[]}, 1, 10); 
    for iunit =  idx_responses'
        for itrial = 1 : length(all_mice_spikes_table{iunit, icondition})
           
            idx = intersect(find(all_mice_spikes_table{iunit, icondition}{itrial}>0), find(all_mice_spikes_table{iunit, icondition}{itrial}<0.5)); 
            X{itrial} = [X{itrial}; all_mice_spikes_table{iunit, icondition}{itrial}(idx)]; 
            %X = [X; all_mice_spikes_table{iunit, icondition}{itrial}(idx)]; 
        end
    end
    for itrial = 1 : length(X)
        if ~isempty(X{itrial})
            %Fit gaussian to this trial
            options = statset('MaxIter',10000);
            gmmModel = fitgmdist(X{itrial}, nGauss, 'Start', 'plus', 'CovarianceType', 'full', 'Options', options);
            pdfValues = pdf(gmmModel, xRange');
            mixingProportions = gmmModel.ComponentProportion;  % Mixing proportions
            means = gmmModel.mu;  % Means of the components
            covariances = gmmModel.Sigma;
           
            [v, idx_ordered_gauss] = sort(means); 
            index_gaussian=0; 
            vect = []; 
            for i =idx_ordered_gauss(1:2)'
                index_gaussian = index_gaussian+1; 
                vect = [vect, means(i), sqrt(covariances(i))]; 
            end
            gfeats = [gfeats; vect([2,4])]; 
            concentration_label = [concentration_label; iodor, concentration]; 
        end
    end
end

size(gfeats)

% Reading and plotting mean intensity curves
params = mean(xlsread('/Users/barrab01/Documents/Repos/ephysMH/perceived_intensity_curves.xlsx')); 

%pure_ppm_odors = {'ethyltiglate':5617.105263, '2-heptanone': 5065.789474, 'ethylbutyrate':16842.10526, 'acetophenone': 526.16}
% ET, 2-HEP, EB, ACE : this is the order of parameters in the perceived
% intensity curves
% odorlabels = {'ethyltiglate', 'ethylbutyrate', 'acetophenone', 'heptanone'}; 
% This is the order in the ephys data labels and colors
parammat = [];
parammat(1,:) = params(1:3); % Ethyltiglate
parammat(2,:) = params(7:9); % Ethyl butyrate
parammat(3,:) = params(10:12); % Acetophenone
parammat(4,:) = params(4:6); % 2-Heptanone
ppm_col = [5617.105263; 16842.10 ; 526.16; 5065.789474]; 
% Colors
int_curves_cols = ['r', 'b', 'y', 'g']; 
odor_symbols = ['o', 'x', 's', 'v']; 

parammat = [parammat, ppm_col]; 
figure
hold on
for io = 1 : 4
    x = logspace(-1, 6, 100); 
    x = log10(x); 
    y = psychocurve(x, parammat(io, 1),  parammat(io, 2), parammat(io, 3)); 
    
    plot(x, y, int_curves_cols(io),  'Linewidth', 2) 
    %set(gca,'Xscale', 'log')
end

%
for ii = 1 : size(concentration_label, 1)
    % convert concentration ot intensity
    io = concentration_label(ii, 1); 
    disp(log10(concentration_label(ii, 2)*parammat(io, 4)))
    realI(ii) = psychocurve(log10(concentration_label(ii, 2)*parammat(io, 4)), parammat(io, 1), parammat(io, 2), parammat(io, 3)); 
    realC(ii) = log10(concentration_label(ii, 2)*parammat(io, 4));  
    plot(log10(concentration_label(ii, 2)*parammat(io, 4)), realI(ii), 'o', 'MarkerFaceColor', int_curves_cols(io)) 

end

% Create data structures and train models
idx_nonnegative = find(concentration_label(:,2)~=0); 

% Intensity prediction
Y = realI(idx_nonnegative)'; 
X = [gfeats(idx_nonnegative,:)]; 
lab_odor = concentration_label(idx_nonnegative,1); 
figure;
hold on; 
MSE_realInt = []; 
for i = 1 : 5
    [rfModel, X_test, Y_test, Y_pred, mse] = emh_random_forest_regressor(X, Y, lab_odor, odor_symbols, 'log'); 
    MSE_realInt(i) = mse; 
end
% Concentration prediction
Y = realC(idx_nonnegative)'; 
X = [gfeats(idx_nonnegative,:)]; 
figure;
hold on; 
MSE_realInt = []; 
for i = 1 : 5
    [rfModel, X_test, Y_test, Y_pred, mse] = emh_random_forest_regressor(X, Y, lab_odor, odor_symbols, 'log'); 
    MSE_realInt(i) = mse; 
end


%% Intensity decoder from gaussians built on firing rates peak latencies 
%  Different mice (same conditions) are separated into different data points
% because we need all trials to compute firing rates
% Read the spikes features
x = [-0.1:0.001:0.5];
xRange = x; 
nGauss = 2; 
gfeats_FR = []; 
conditions = [2:8, 10:16, 18:24, 26:32]; % skipping no stim
concentration_label = []; 
for icondition = conditions
    % Find parameters
    iodor = all_mice_condition_table(icondition, 1); 
    iconc = all_mice_condition_table(icondition, 2); 
    concentration = all_mice_condition_table(icondition, 3); 
    
    % Filter units that have a peak
    idx_responses  =find(~isnan(peak_latency(:, icondition))); 
    
    % Extract all spikes
    for im = 1 : length(unique(imouse)) % separate different mice
        idx_mouse = find(imouse == im); 
        idx_units_thismouse = intersect(idx_responses, idx_mouse); 
        
        X = peak_latency(idx_units_thismouse, icondition); 
        %X = [X; all_mice_spikes_table{iunit, icondition}{itrial}(idx)]; 
        options = statset('MaxIter',10000);
        gmmModel = fitgmdist(X, nGauss, 'Start', 'plus', 'CovarianceType', 'full', 'Options', options);
        pdfValues = pdf(gmmModel, xRange');
        mixingProportions = gmmModel.ComponentProportion;  % Mixing proportions
        means = gmmModel.mu;  % Means of the components
        covariances = gmmModel.Sigma;
       
        [v, idx_ordered_gauss] = sort(means); 
        index_gaussian=0; 
        vect = []; 
        for i =idx_ordered_gauss(1:2)' % only the first 2 gaussians
            index_gaussian = index_gaussian+1; 
            vect = [vect, means(i), sqrt(covariances(i))]; 
        end
        gfeats_FR = [gfeats_FR; vect]; 
        concentration_label = [concentration_label; iodor, concentration]; 
    end
end


size(gfeats_FR)

% Reading and plotting mean intensity curves
params = mean(xlsread('/Users/barrab01/Documents/Repos/ephysMH/perceived_intensity_curves.xlsx')); 

%pure_ppm_odors = {'ethyltiglate':5617.105263, '2-heptanone': 5065.789474, 'ethylbutyrate':16842.10526, 'acetophenone': 526.16}
% ET, 2-HEP, EB, ACE : this is the order of parameters in the perceived
% intensity curves
% odorlabels = {'ethyltiglate', 'ethylbutyrate', 'acetophenone', 'heptanone'}; 
% This is the order in the ephys data labels and colors
parammat = [];
parammat(1,:) = params(1:3); % Ethyltiglate
parammat(2,:) = params(7:9); % Ethyl butyrate
parammat(3,:) = params(10:12); % Acetophenone
parammat(4,:) = params(4:6); % 2-Heptanone
ppm_col = [5617.105263; 16842.10 ; 526.16; 5065.789474]; 
% Colors
int_curves_cols = ['r', 'b', 'y', 'g']; 
odor_symbols = ['o', 'x', 's', 'v']; 

parammat = [parammat, ppm_col]; 
figure
hold on
for io = 1 : 4
    x = logspace(-1, 6, 100); 
    x = log10(x); 
    y = psychocurve(x, parammat(io, 1),  parammat(io, 2), parammat(io, 3)); 
    
    plot(x, y, int_curves_cols(io),  'Linewidth', 2) 
    %set(gca,'Xscale', 'log')
end

%
realI = [];
realC = [];
svpC = [];
for ii = 1 : size(concentration_label, 1)
    % convert concentration ot intensity
    io = concentration_label(ii, 1); 
    disp(log10(concentration_label(ii, 2)*parammat(io, 4)))
    realI(ii) = psychocurve(log10(concentration_label(ii, 2)*parammat(io, 4)), parammat(io, 1), parammat(io, 2), parammat(io, 3)); 
    realC(ii) = log10(concentration_label(ii, 2)*parammat(io, 4));  
    svpC(ii) = log10(concentration_label(ii, 2));  
    plot(log10(concentration_label(ii, 2)*parammat(io, 4)), realI(ii), 'o', 'MarkerFaceColor', int_curves_cols(io)) 

end

% Create data structures and train models
idx_nonnegative = find(concentration_label(:,2)~=0); 

% Intensity prediction
Y = realI(idx_nonnegative)'; 
X = [gfeats_FR(idx_nonnegative,:)]; 
lab_odor = concentration_label(idx_nonnegative,1); 
figure;
hold on; 
MSE_realInt = []; 
for i = 1 : 10
    [rfModel, X_test, Y_test, Y_pred, mse] = emh_random_forest_regressor(X, Y, lab_odor, odor_symbols, 'log'); 
    MSE_realInt(i) = mse; 
end
% ppm Concentration prediction
Y = realC(idx_nonnegative)'; 
X = [gfeats_FR(idx_nonnegative,:)]; 
figure;
hold on; 
MSE_ppmconc = []; 
for i = 1 : 10
    [rfModel, X_test, Y_test, Y_pred, mse] = emh_random_forest_regressor(X, Y, lab_odor, odor_symbols, 'log'); 
    MSE_ppmconc(i) = mse; 
end

% SVP Concentration prediction
Y = svpC(idx_nonnegative)'; 
X = [gfeats_FR(idx_nonnegative,:)]; 
figure;
hold on; 
MSE_svpconc = []; 
for i = 1 : 10
    [rfModel, X_test, Y_test, Y_pred, mse] = emh_random_forest_regressor(X, Y, lab_odor, odor_symbols, 'log'); 
    MSE_svpconc(i) = mse; 
end

figure
boxplot([MSE_realInt, MSE_ppmconc, MSE_svpconc], [1*ones(size(MSE_realInt)), 2*ones(size(MSE_ppmconc)), 3*ones(size(MSE_svpconc))])
ranksum(MSE_ppmconc, MSE_realInt)

%% Neural network based model

% Example Data (regression task: predicting y from x)
%T = Y
% SVP Concentration prediction
scale = 'log'; 
Y = realI(idx_nonnegative)'; 
X = [gfeats_FR(idx_nonnegative,:)]; 
figure
% Create a feedforward neural network (1 hidden layer with 10 neurons)
NNeurons = 10; 
trainRatio = 0.8;
transferFun = 'purelin'; 
for i = 1 : 10
    [net, X_test, Y_test, Y_pred, nrmse] = emh_perceptron_regressor(X, Y, lab_odor, odor_symbols, scale, NNeurons, trainRatio, transferFun); 
    MSE_realInt(i) = nrmse; 
end

Y = realC(idx_nonnegative)'; 
X = [gfeats_FR(idx_nonnegative,:)]; 
figure
% Create a feedforward neural network (1 hidden layer with 10 neurons)
NNeurons = 10; 
trainRatio = 0.8;
transferFun = 'purelin'; 
for i = 1 : 10
    [net, X_test, Y_test, Y_pred, nrmse] = emh_perceptron_regressor(X, Y, lab_odor, odor_symbols, scale, NNeurons, trainRatio, transferFun); 
    MSE_ppmconc(i) = nrmse; 
end


Y = svpC(idx_nonnegative)'; 
X = [gfeats_FR(idx_nonnegative,:)]; 
figure
% Create a feedforward neural network (1 hidden layer with 10 neurons)
NNeurons = 10; 
trainRatio = 0.8;
transferFun = 'purelin'; 
for i = 1 : 10
    [net, X_test, Y_test, Y_pred, nrmse] = emh_perceptron_regressor(X, Y, lab_odor, odor_symbols, scale, NNeurons, trainRatio, transferFun); 
    MSE_svpconc(i) = nrmse; 
end

figure
boxplot([MSE_realInt, MSE_ppmconc, MSE_svpconc], [1*ones(size(MSE_realInt)), 2*ones(size(MSE_ppmconc)), 3*ones(size(MSE_svpconc))])
ranksum(MSE_ppmconc, MSE_realInt)


%% Testing peak finding
%[FRmatthismouse, peak_latencythismouse] = extract_peak_latency(spikes_table, sig, mint, maxt); 

%% Checking peaks of firing rate curves

figureC= figure()
ms100 = odor_onset + abs(0.1/sig); 
ms200 = odor_onset + abs(0.2/sig); 
% Sort and plot all mice together 
% sort responses across all mice
iconc = 5; 
for iodor = 1 : Nodors % for each odor
    icondition = intersect(find(condition_table(:,1) == iodor), find(condition_table(:,2) == iconc)); 
    % sort for each condition separately
    [sortedpeaks, sorted_idx{iodor}] = sort(peak_latency(:, icondition)); 
end
for icondition = 1 : size(peak_latency, 2)
    
    
    iodor = condition_table(icondition, 1); 
    iconc = condition_table(icondition, 2); 
    % sort for each condition separately
    [sortedpeaks, sorted_idx{iodor}] = sort(peak_latency(:, icondition));  % each
    %condition per itself

    ii = 0; 
    for i = sorted_idx{iodor}'
        ii = ii + 1; 
        sorted_responses(ii, :) = FRmat{icondition}(i, :); 
        newpeaks(ii) = peak_latency(i, icondition); 
    end
    
    % Firing rates curves
    figure(figureC)
    subplot(Nodors+2, Nconcs, (iodor-1)*Nconcs + iconc)
    tresp = linspace(0, sig*size(sorted_responses, 2), size(sorted_responses, 2)); 
    for iresp = 1 : size(sorted_responses, 1)
        plot(tresp, sorted_responses(iresp,:))
        hold on
        plot(newpeaks(iresp), max(sorted_responses(iresp,:)), 'o')
        pause
        clf
    end
end
%% Compute psths
sig = 0.010; 
mint = 0; 
maxt = 0.5; 
tvect = linspace(mint, maxt, (maxt-mint)/sig); 
figure
for icondition = 1 : size(spikes_table, 2)
   
    FRmat = []; 
    sorted_responses = []; 
    for iunit = 1 : size(spikes_table, 1)
        [FR,t, err] = compute_FR_from_raster(spikes_table{iunit, icondition} , sig, mint, maxt, tvect ); 
        FRmat = [FeRmat; FR]; 
        [xmax] = findpeaks(FR, 5); 
        %[maxv, maxloc] = max(FR(xmax.loc)); 
        [maxidx]= find_true_peak(FR, xmax); 
        peak_latency(iunit) = maxidx*sig; 
        
%         if ~isnan(maxidx)
%             plot(t, FR)
%             hold on
%             plot(t(maxidx), FR(maxidx), 'or')
%             pause
%         end
    end
    [sortedpeaks, sorted_idx] = sort(peak_latency); 
    ii = 0; 
    for i = sorted_idx
        ii = ii + 1; 
        sorted_responses(ii, :) = FRmat(i, :); 
    end
    
    iodor = condition_table(icondition, 1); 
    iconc = condition_table(icondition, 2); 
    subplot(Nodors, Nconcs, (iodor-1)*Nconcs + iconc)
    imagesc((sorted_responses), [0,50])
    hold on
    idx_inh = floor(size(sorted_responses, 2)*abs(-mint/(-mint + maxt))); 
    plot([idx_inh, idx_inh], [0, size(sorted_responses, 1)], '--w', 'Linewidth', 2) 
end


%% SINGLE MOUSE
%% Extract spikes from each trials
Nodors = 4; 
Nconcs = 8; 
mywindow = 0.5; 
spikes_table = cell(length(PC_units_idx), Nodors*Nconcs);
count_table = zeros(length(PC_units_idx), Nodors*Nconcs); 
condition_table = zeros( Nodors*Nconcs, 2); 
for itrial  = 1 : length(Odors)
    
    % Find odor and concentration indexes to pull together trials with the
    % same conditions
    iodor = find(all_odors == Odors(itrial)); 
    idx_odor = find(Odors ==iodor+4); 
    all_conc_thisodor = unique(Concentration(idx_odor));
    iconc = find(all_conc_thisodor == Concentration(itrial)); 
    condition_idx = (iodor-1)*length(all_conc_thisodor) + iconc; 
    
    condition_table( condition_idx, :) = [iodor, iconc]; 

    all_prex_before_fvo_idx = find(Prex<FVO(itrial));
    all_prex_after_fvo_idx = find(Prex>FVO(itrial));
    all_postx_before_fvo_idx=  find(Postx<FVO(itrial));
    
    % Start and end of baseline
    start_bas = Prex(all_prex_before_fvo_idx(end-2));  
    end_bas = start_bas + mywindow; Postx(all_postx_before_fvo_idx(end)); 
    
    % Start and end of inh
    start_inh = Prex(all_prex_after_fvo_idx(1)); 
    end_inh =start_inh + mywindow; 
    
    % Loop through units and find spikes in baseline and first sniff
    iiunit = 0; 
    for iunit = PC_units_idx'
        iiunit = iiunit+1; 
        last_spike = nwbFile.units.spike_times_index.data(iunit); 
        if iunit >1
            first_spike = nwbFile.units.spike_times_index.data(iunit-1)+1; 
        else
            first_spike = 1; 
        end

        allspikes = spikes(first_spike:last_spike); 
        if isempty(spikes_table{iiunit, condition_idx})
            spikes_table{iiunit, condition_idx} = {allspikes-start_inh}; 
            count_table(iiunit, condition_idx) = count_table(iiunit, condition_idx) +1; 
        else
            spikes_table{iiunit, condition_idx}{count_table(iiunit, condition_idx) +1} = allspikes-start_inh; 
            count_table(iiunit, condition_idx) = count_table(iiunit, condition_idx) +1;
        end
    end
    percent_complete = (itrial / length(Odors)) * 100;
    disp([' ----   ', num2str(percent_complete), ' % of Extracting spikes   ----'])
end

%% Compute psths
sig = 0.010; 
mint = 0; 
maxt = 0.5; 
tvect = linspace(mint, maxt, (maxt-mint)/sig); 
figure
for icondition = 1 : size(spikes_table, 2)
   
    FRmat = []; 
    sorted_responses = []; 
    for iunit = 1 : size(spikes_table, 1)
        [FR,t, err] = compute_FR_from_raster(spikes_table{iunit, icondition} , sig, mint, maxt, tvect ); 
        FRmat = [FeRmat; FR]; 
        [xmax] = findpeaks(FR, 5); 
        %[maxv, maxloc] = max(FR(xmax.loc)); 
        [maxidx]= find_true_peak(FR, xmax); 
        peak_latency(iunit) = maxidx*sig; 
        
%         if ~isnan(maxidx)
%             plot(t, FR)
%             hold on
%             plot(t(maxidx), FR(maxidx), 'or')
%             pause
%         end
    end
    [sortedpeaks, sorted_idx] = sort(peak_latency); 
    ii = 0; 
    for i = sorted_idx
        ii = ii + 1; 
        sorted_responses(ii, :) = FRmat(i, :); 
    end
    
    iodor = condition_table(icondition, 1); 
    iconc = condition_table(icondition, 2); 
    subplot(Nodors, Nconcs, (iodor-1)*Nconcs + iconc)
    imagesc((sorted_responses), [0,50])
    hold on
    idx_inh = floor(size(sorted_responses, 2)*abs(-mint/(-mint + maxt))); 
    plot([idx_inh, idx_inh], [0, size(sorted_responses, 1)], '--w', 'Linewidth', 2) 
end
