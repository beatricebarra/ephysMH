%% Script to produce figures for Achems poster 2025

%% Plot units repartition for each mouse
savefilenameallmice = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/all_mice_dataset_OB.mat'; 
load(savefilenameallmice)
figure
for mouse_i= 1:7
    all_units = find(imouse == mouse_i); 
    
    % Find which units are responsive for each odor in this mouse
    [selunits] = select_responsive_units(all_mice_spikes_table, all_mice_condition_table, all_units); 
    
    % Find which units are responsive at least once
    sometimes_responsive = selunits{1}; 
    for iodor = 2 : 4
        sometimes_responsive =union(sometimes_responsive, selunits{iodor});  
    end
    % Find which units are never responsive
    never_responsive = length(all_units) - length(sometimes_responsive); 
    % Find units that are responsive to multiple odors 
    multiple_responses = length(all_units) - never_responsive;
    % Find exclusive responses for each odorant
    for iodor = 1 : 4
        exclusive_responses{iodor} = selunits{iodor}; 
        other_odors =setdiff([1:1:4], [iodor]); 
        for oodor = other_odors
            exclusive_responses{iodor} =setdiff(exclusive_responses{iodor}, selunits{oodor});  
        end
        multiple_responses = multiple_responses-length(exclusive_responses{iodor}); 
    end
    
    
    subplot(1, 7, mouse_i)

    pie([ never_responsive, multiple_responses,...
        length(exclusive_responses{1}), length(exclusive_responses{2}), ...
        length(exclusive_responses{3}), length(exclusive_responses{4})])
    pielabels = {'No resp', 'Many', 'odor1', 'odor2', 'odor3', 'odor4'}
    legend(pielabels)
    sum([multiple_responses, never_responsive, length(exclusive_responses{1}), length(exclusive_responses{2}), length(exclusive_responses{3}), length(exclusive_responses{4})])

end

%% Plot firing rates heatmaps

clear sorted_responses
Nodors = 4
Nconcs = 8; 
odorcols = ['r', 'b', 'y', 'g']; 
dt_FR = 0.010; 
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
clear sorted_responses
mouse_label =  []; 
odor_label =  []; 
concentration = []; 
ppm_concentration =  []; 
intensity = []; 
histdata = []; 
frdata = []; 
for mouse_i= 1:7
    all_units = find(imouse == mouse_i); 
    figureFR = figure()
    figureH= figure()
    figtraces = figure()
    actual_sorted_units= {}; 
    % Find which units are responsive for each odor in this mouse
    [selunits] = select_responsive_units(all_mice_spikes_table, all_mice_condition_table, all_units); 
    iconditions = (mouse_i-1)*32+1:(mouse_i)*32;
    % Sort at a specific condition for each odor
    NConcs = length(unique(unique(all_mice_condition_table(:, 2)))); 
    for iodor = unique(all_mice_condition_table(:, 1))'
        iconc = 5; 
        icondition = (iodor-1)*NConcs + iconc; 
        [sortedpeaks, sorted_idx{iodor}] = sort(peak_latency(selunits{iodor}, icondition));  
        actual_sorted_units{iodor} = selunits{iodor}(sorted_idx{iodor}); 
    end
    for icondition = 1: length(iconditions)
        % incondition is the index across all conditions in the
        % "all_mice_condition_table" while icond is an index ranging 1:32
        % that specifies the actual condition which is the same for all
        % mice 
        icond = iconditions(icondition); 
        % Find odor and conceentration
        iodor = all_mice_condition_table(icond, 1); 
        iconc = all_mice_condition_table(icond, 2); 
        conc = all_mice_condition_table(icond, 3); 
        nreps = all_mice_condition_table(icond, 4); 
        
        % Sort for each condition separately and plot the units subsequence
        % for only the units implicated in this odor in this mouse
        %[sortedpeaks, sorted_idx{iodor}] = sort(peak_latency(selunits{iodor}, icondition));  
        %actual_sorted_units = selunits{iodor}(sorted_idx{iodor}); 
        ii = 0; 
        newpeaks = []; 
        sorted_responses = []; 
        for i = actual_sorted_units{iodor}
            ii = ii + 1; 
            if ~isnan(peak_latency(i, icondition)) %Check that there is peak
                sorted_responses(ii, :) = allmice_FRmat{icondition}(i, :); 
                newpeaks(ii) = peak_latency(i, icondition);
                FRpeaks(ii) = sorted_responses(ii, round(newpeaks(ii)/dt_FR)+1); 
                avgFR(ii) = mean(allmice_FRmat{icondition}(i, 1:20)); 
            else
                FRpeaks(ii) = 0; 
                avgFR(ii) = mean(allmice_FRmat{icondition}(i, 1:20)); 
            end
        end
        
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
        bin_centers = linspace(0.01, 0.5, 50); 
        hist(newpeaks, bin_centers);
        [trialhistdata] = hist(newpeaks, bin_centers); 
        title([num2str(iodor) , num2str(iconc) , ])
        ylim([0, 80])
        
        % Building dataset for regressor on histogram
        realI = computeIntensity(intensity_file, [iodor, conc]); 
        mouse_label = [mouse_label; mouse_i]; 
        odor_label = [odor_label; iodor]; 
        concentration = [concentration; conc]; 
        ppm_concentration = [ppm_concentration; conc*ppm_col(iodor)]; 
        intensity = [intensity; realI]; 
        histdata = [histdata; trialhistdata]; 
        frdata = [frdata; avgFR(1:100)]; 
        % 
        figure(figtraces)
        subplot(1, 4, iodor)
        hold on
        if iconc ==1
             plot(mean((sorted_responses)), 'Color', 'k')
        else
            plot(mean((sorted_responses)), 'Color', yellow2redmap(iconc,:))
        end
        ylim([0, 30])
       

    end
    
%     figure(figureFR)
%     figname = ['/Users/barrab01/Dropbox/CONFERENCES/2025/Achems/Poster/FRheatmaps_M', num2str(mouse_i), '.pdf']; 
%     saveas(gcf,figname)
%     
%     figure(figureH)
%     figname = ['/Users/barrab01/Dropbox/CONFERENCES/2025/Achems/Poster/histograms_M', num2str(mouse_i), '.pdf']; 
%     saveas(gcf,figname)
%     
%     figure(figtraces)
%     figname = ['/Users/barrab01/Dropbox/CONFERENCES/2025/Achems/Poster/Compound_reponse_M', num2str(mouse_i), '.pdf']; 
%     saveas(gcf,figname)

end

%% Find time window for better accuracy from time data
windows = [(1:5:46)', (5:5:50)']; 
for iw = 1 : size(windows)
    
    ws = windows(iw, 1); 
    wss = windows(iw, 2); 
    bincenter(iw) = mean([bin_centers(ws), bin_centers(wss)]); 
    idx_nonnegative = (find(concentration~=0)); %Avoid concentration = 0 because logarithm becomes inf
    Y = [log10(intensity(idx_nonnegative)), (log10(concentration(idx_nonnegative))),  (log10(ppm_concentration(idx_nonnegative)))]; 
    X = [histdata(idx_nonnegative, ws:wss)]; 

    X_input = X;      % 2 x N
    Y_target = Y(:, 1);     % 1 x N   
    xfit = linspace(min(Y_target), max(Y_target), 100); 
    for i =  1 : 50
        % SVM
        % Create testing and training set for each run
        trainRatio = 0.8;
        idx = randperm(length(Y_target));
        trainIdx = idx(1:round(trainRatio * length(Y_target)));
        testIdx = idx(round(trainRatio * length(Y_target)) + 1:end);
        
        X_train = X_input(trainIdx, :);
        Y_train = Y_target(trainIdx);
        X_test = X_input(testIdx, :);
        Y_test = Y_target(testIdx);
         
        svmModel = fitrsvm(X_train, Y_train, ...
            'KernelFunction', 'gaussian', ...
            'KernelScale', 'auto', ...
            'Standardize', true);
        Y_pred = predict(svmModel, X_test);
        
        % Plot actual vs predicted
        plot(Y_test, Y_pred, 'o');
        xlabel('True Y'); ylabel('Predicted Y');
        title('Regression: True vs Predicted');
        grid on; axis equal;
        hold on
        plot(Y_test, Y_test, 'r', 'Linewidth', 2)
        
        linearCoef = polyfit(Y_test,Y_pred,1);
        linearFit(i,:) = polyval(linearCoef,xfit)';
        
        % Compute mean squared error
        mse_error = mse(Y_pred - Y_test);
        rmse = sqrt(median((Y_test - Y_pred).^2));
        nrmse = rmse/(max([Y_test]) - min([Y_test]));
        net_errors_window(i , iw) = nrmse; 
        
        % Comput ePearsons coefficient
        R = corrcoef(Y_pred, Y_test);
        Rcoeff_window(i , iw)  = R(1,2);  % This is the Pearson correlation coefficient
    end
end

figure
errorbar(bincenter, mean(net_errors_window), std(net_errors_window), 'o-')
hold on
errorbar(bincenter, mean(Rcoeff_window), std(Rcoeff_window), 'o-')

%% Run regressor model on histogram data

idx_nonnegative = (find(concentration~=0)); %Avoid concentration = 0 because logarithm becomes inf
Y = [log10(intensity(idx_nonnegative)), (log10(concentration(idx_nonnegative))),  (log10(ppm_concentration(idx_nonnegative)))]; 
X = [histdata(idx_nonnegative, 5:15)]; 

for ioutput = 1 : size(Y, 2) % Those are different quantities to predict 
    figure;
    % SVM
    X_input = X;      % 2 x N
    Y_target = Y(:, ioutput);     % 1 x N   
    xfit = linspace(min(Y_target), max(Y_target), 100); 
    for i =  1 : 50
        % SVM
        
        % Create testing and training set for each run
        trainRatio = 0.8;
        idx = randperm(length(Y_target));
        trainIdx = idx(1:round(trainRatio * length(Y_target)));
        testIdx = idx(round(trainRatio * length(Y_target)) + 1:end);
        
        X_train = X_input(trainIdx, :);
        Y_train = Y_target(trainIdx);
        X_test = X_input(testIdx, :);
        Y_test = Y_target(testIdx);
         
        svmModel = fitrsvm(X_train, Y_train, ...
            'KernelFunction', 'gaussian', ...
            'KernelScale', 'auto', ...
            'Standardize', true);
        Y_pred = predict(svmModel, X_test);
        
        % Plot actual vs predicted
        plot(Y_test, Y_pred, 'o');
        xlabel('True Y'); ylabel('Predicted Y');
        title('Regression: True vs Predicted');
        grid on; axis equal;
        hold on
        plot(Y_test, Y_test, 'r', 'Linewidth', 2)
        
        linearCoef = polyfit(Y_test,Y_pred,1);
        linearFit(i,:) = polyval(linearCoef,xfit)';
        
        % Compute mean squared error
        mse_error = mse(Y_pred - Y_test);
        rmse = sqrt(median((Y_test - Y_pred).^2));
        nrmse = rmse/(max([Y_test]) - min([Y_test]));
        net_errors(i , ioutput) = nrmse; 
        
        % Comput ePearsons coefficient
        R = corrcoef(Y_pred, Y_test);
        Rcoeff(i , ioutput)  = R(1,2);  % This is the Pearson correlation coefficient
    end
    hold on
    %plot(Y_test,linearFit,'k-', 'Linewidth', 2)
    
    stdshade(linearFit, 0.3, 'k', xfit)
end

figure
subplot(1, 2, 1)
boxplot(net_errors)
set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'Intensity', 'C(svp)', 'C(ppm)'})
title('Normalized RMSE error')
%ylim([0, 0.12])
p1= ranksum(net_errors(:, 1), net_errors(:, 2))
p2= ranksum(net_errors(:, 1), net_errors(:, 3))

subplot(1, 2, 2)
boxplot(Rcoeff)
title('Pearsons coefficient')
set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'Intensity', 'C(svp)', 'C(ppm)'})
ylim([0.5, 1])
p1= ranksum(Rcoeff(:, 1), Rcoeff(:, 2))
p2= ranksum(Rcoeff(:, 1), Rcoeff(:, 3))

%% Run regressor model on firing rate data

idx_nonnegative = (find(concentration~=0)); %Avoid concentration = 0 because logarithm becomes inf
Y = [log10(intensity(idx_nonnegative)), (log10(concentration(idx_nonnegative))),  (log10(ppm_concentration(idx_nonnegative)))]; 
X = [frdata(idx_nonnegative,:)]; 

for ioutput = 1 : size(Y, 2) % Those are different quantities to predict 
    figure;
    % SVM
    X_input = X;      % 2 x N
    Y_target = Y(:, ioutput);     % 1 x N   
    xfit = linspace(min(Y_target), max(Y_target), 100); 
    for i =  1 : 50
        % SVM
        
        % Create testing and training set for each run
        trainRatio = 0.8;
        idx = randperm(length(Y_target));
        trainIdx = idx(1:round(trainRatio * length(Y_target)));
        testIdx = idx(round(trainRatio * length(Y_target)) + 1:end);
        
        X_train = X_input(trainIdx, :);
        Y_train = Y_target(trainIdx);
        X_test = X_input(testIdx, :);
        Y_test = Y_target(testIdx);
         
        svmModel = fitrsvm(X_train, Y_train, ...
            'KernelFunction', 'gaussian', ...
            'KernelScale', 'auto', ...
            'Standardize', true);
        Y_pred = predict(svmModel, X_test);
        
        % Plot actual vs predicted
        plot(Y_test, Y_pred, 'o');
        xlabel('True Y'); ylabel('Predicted Y');
        title('Regression: True vs Predicted');
        grid on; axis equal;
        hold on
        plot(Y_test, Y_test, 'r', 'Linewidth', 2)
        
        linearCoef = polyfit(Y_test,Y_pred,1);
        linearFit(i,:) = polyval(linearCoef,xfit)';
        
        % Compute mean squared error
        mse_error = mse(Y_pred - Y_test);
        rmse = sqrt(median((Y_test - Y_pred).^2));
        nrmse = rmse/(max([Y_test]) - min([Y_test]));
        net_errors(i , ioutput) = nrmse; 
        
        % Comput ePearsons coefficient
        R = corrcoef(Y_pred, Y_test);
        Rcoeff(i , ioutput)  = R(1,2);  % This is the Pearson correlation coefficient
    end
    hold on
    %plot(Y_test,linearFit,'k-', 'Linewidth', 2)
    
    stdshade(linearFit, 0.3, 'k', xfit)
end

figure
subplot(1, 2, 1)
boxplot(net_errors)
set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'Intensity', 'C(svp)', 'C(ppm)'})
title('Normalized RMSE error')
%ylim([0, 0.12])
p1= ranksum(net_errors(:, 1), net_errors(:, 2))
p2= ranksum(net_errors(:, 1), net_errors(:, 3))

subplot(1, 2, 2)
boxplot(Rcoeff)
title('Pearsons coefficient')
set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'Intensity', 'C(svp)', 'C(ppm)'})
ylim([0.5, 1])
p1= ranksum(Rcoeff(:, 1), Rcoeff(:, 2))
p2= ranksum(Rcoeff(:, 1), Rcoeff(:, 3))
%% Plot example psth for some units


%% Read complete synchrony dataset with all mice
synch_file = '/Users/barrab01/Documents/Repos/ephysMH/processed_data/synchrony_timeseries_OB.mat'; 
synch_condition_file = '/Users/barrab01/Documents/Repos/ephysMH/processed_data/synch_condition_OB.mat'; 

load(synch_file)
load(synch_condition_file)

mintime = -0.2; 
maxtime = 0.5; 
dt = 0.001; 
maxW = 0.01; 
%% Plot synchrony measure profiles, averaged across condition, for all mice, all odors, all concentrations
mycolormaps = generateColorMaps(10); 
colorlabel = {'red', 'blue', 'yellow', 'green'}; 
figure
for mouse_i = 1: 7
    theodors = unique(synch_condition{mouse_i}(:, 2));

    for iodor = 1 : length(unique(theodors))
        idxO = find(synch_condition{mouse_i}(:,2) == theodors(iodor)); 
        concentrations = unique(synch_condition{mouse_i}(idxO,1)); 
        for ic = [3,5,7]%1 : length(concentrations)
            idxC = find(synch_condition{mouse_i}(:,1) == concentrations(ic)); 
            idx = intersect(idxO, idxC); 
            ax(mouse_i, iodor) = subplot(7, length(unique(theodors)), (mouse_i-1)*4 + iodor); 
            hold on
            if length(idx)> 1
                stdshade(synchrony_timeseries{mouse_i}(idx, :), 0.3, mycolormaps.(colorlabel{iodor})(ic,:))
            end
            avg_profile = mean(synchrony_timeseries{mouse_i}(idx, :)); 
        end
    end
end

linkaxes(ax, 'x')

%% Create dataset for regressor model that includes peak value and location of synchrony time_series as computed with my measure
% Repetitions are averaged to create a single datapoint
odor_symbols = ['o', 's', 't', 'x']; 
intensity_file  = '/Users/barrab01/Documents/Repos/ephysMH/perceived_intensity_curves.xlsx'; 

time_points = [mintime:dt:maxtime]
ppm_col = [5617.105263; 16842.10 ; 526.16; 5065.789474];
X = []; 
Y = []; 
mouse_label = []; 
odor_label = []; 
intensity = []; 
concentration = []; 
synch = []; 
ppm_concentration = []; 
time_idx_ofinterest = (intersect(find(time_points>=0), find(time_points<=0.2))); 
theodors = unique(synch_condition{mouse_i}(:,2));
for mouse_i = 1: 7
    for iodor = [1, 4] %1 : length(unique(theodors))
        idxO = find(synch_condition{mouse_i}(:,2) == theodors(iodor)); 
        concentrations = unique(synch_condition{mouse_i}(idxO,1)); 
        for ic = 1 : length(concentrations)
            idxC = find(synch_condition{mouse_i}(:,1) == concentrations(ic)); 
            idx = intersect(idxO, idxC); 
            realI = computeIntensity(intensity_file, [iodor, concentrations(ic)]); 
        
            mouse_label = [mouse_label; mouse_i]; 
            odor_label = [odor_label; iodor]; 
            concentration = [concentration; concentrations(ic)]; 
            ppm_concentration = [ppm_concentration; concentrations(ic)*ppm_col(iodor)]; 
            intensity = [intensity; realI]; 
            
            avg_profile = mean(synchrony_timeseries{mouse_i}(idx, :), 1); 
            [maxv, maxloc] = max(avg_profile(time_idx_ofinterest)');
            synch = [synch;  maxv, time_points(time_idx_ofinterest(maxloc))]; %maxv',%time_points(time_idx_ofinterest(maxloc))'
        end
    end
end


%  Apply svm regressor model 
close all

idx_nonnegative = (find(concentration~=0)); %Avoid concentration = 0 because logarithm becomes inf
Y = [log10(intensity(idx_nonnegative)), (log10(concentration(idx_nonnegative))),  (log10(ppm_concentration(idx_nonnegative)))]; 
X = [synch(idx_nonnegative, :)]; 

for ioutput = 1 : size(Y, 2) % Those are different quantities to predict 
    figure;
    % SVM
    X_input = X;      % 2 x N
    Y_target = Y(:, ioutput);     % 1 x N   
    xfit = linspace(min(Y_target), max(Y_target), 100); 
    for i =  1 : 50
        % SVM
        
        % Create testing and training set for each run
        trainRatio = 0.8;
        idx = randperm(length(Y_target));
        trainIdx = idx(1:round(trainRatio * length(Y_target)));
        testIdx = idx(round(trainRatio * length(Y_target)) + 1:end);
        
        X_train = X_input(trainIdx, :);
        Y_train = Y_target(trainIdx);
        X_test = X_input(testIdx, :);
        Y_test = Y_target(testIdx);
         
        svmModel = fitrsvm(X_train, Y_train, ...
            'KernelFunction', 'linear', ...
            'KernelScale', 'auto', ...
            'Standardize', true);
        Y_pred = predict(svmModel, X_test);
        
        % Plot actual vs predicted
        plot(Y_test, Y_pred, 'o');
        xlabel('True Y'); ylabel('Predicted Y');
        title('Regression: True vs Predicted');
        grid on; axis equal;
        hold on
        plot(Y_test, Y_test, 'r', 'Linewidth', 2)
        
        linearCoef = polyfit(Y_test,Y_pred,1);
        linearFit(i,:) = polyval(linearCoef,xfit)';
        
        % Compute mean squared error
        mse_error = mse(Y_pred - Y_test);
        rmse = sqrt(median((Y_test - Y_pred).^2));
        nrmse = rmse/(max([Y_test]) - min([Y_test]));
        net_errors(i , ioutput) = nrmse; 
        
        % Comput ePearsons coefficient
        R = corrcoef(Y_pred, Y_test);
        Rcoeff(i , ioutput)  = R(1,2);  % This is the Pearson correlation coefficient
    end
    hold on
    %plot(Y_test,linearFit,'k-', 'Linewidth', 2)
    
    stdshade(linearFit, 0.3, 'k', xfit)
end

figure
subplot(1, 2, 1)
boxplot(net_errors)
set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'Intensity', 'C(svp)', 'C(ppm)'})
title('Normalized RMSE error')
%ylim([0, 0.12])
p1= ranksum(net_errors(:, 1), net_errors(:, 2))
p2= ranksum(net_errors(:, 1), net_errors(:, 3))

subplot(1, 2, 2)
boxplot(Rcoeff)
title('Pearsons coefficient')
set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'Intensity', 'C(svp)', 'C(ppm)'})
ylim([0.5, 1])
p1= ranksum(Rcoeff(:, 1), Rcoeff(:, 2))
p2= ranksum(Rcoeff(:, 1), Rcoeff(:, 3))

%% Create dataset for regressor model that includes peak value and location of synchrony time_series as computed with my measure
% Repetitions are averaged to create a single datapoint

% Repetitions are averaged to create a single datapoint
odor_symbols = ['o', 's', 't', 'x']; 
intensity_file  = '/Users/barrab01/Documents/Repos/ephysMH/perceived_intensity_curves.xlsx'; 
ppm_col = [5617.105263; 16842.10 ; 526.16; 5065.789474];
X = []; 
Y = []; 
intensity_class = []; 
synch_class = {[], []}; 


ppm_concentration_class = []; 
ppm_concentration_values = [];
time_idx_ofinterest = (intersect(find(time_points>=0), find(time_points<=0.2))); 
temp = []; 
for mouse_i = 1: 7
    
    for iodor = 1 : length(unique(theodors))
        idxO = find(synch_condition{mouse_i}(:,2) == theodors(iodor)); 
        concentrations = unique(synch_condition{mouse_i}(idxO,1)); 
        
        for ic = 3 : 6 % defin concentration classes 
            temp = [temp; log10(concentrations(ic)*ppm_col(iodor))]; 
        end
    end
   
    
end  
step = (1*max(temp) - 1*min(temp)) / 3;
t1 = min(temp) + step;      % Lower threshold
t2 = min(temp) + 2 * step; 
figure, hold on, plot(temp, rand(size(temp)), 'o'), plot([t1, t1], [0, 1]), plot([t2, t2], [0, 1])

plot([min(temp), min(temp)], [0, 1])

odor_picks_ppm = {[3 : 6], [3 : 6], [3 : 5], [3 : 6]}; 
odor_picks = {[4 : 6], [4 : 6], [4 : 5], [4 : 6]}; 

for mouse_i = 1: 7    % Upper threshold
    for iodor = 1 : length(unique(theodors)) % [1, 4]
        idxO = find(synch_condition{mouse_i}(:,2) == theodors(iodor)); 
        concentrations = unique(synch_condition{mouse_i}(idxO,1)); 
        for ic = odor_picks_ppm{iodor}
            idxC = find(synch_condition{mouse_i}(:,1) == concentrations(ic)); 
            idx = intersect(idxO, idxC); 
            ppm_concentration_values = [ppm_concentration_values; log10(concentrations(ic)*ppm_col(iodor))]; 
            
            avg_profile = mean(synchrony_timeseries{mouse_i}(idx, :), 1); 
            [maxv, maxloc] = max(avg_profile(time_idx_ofinterest)');
            synch_class{2} = [synch_class{2};  maxv, time_points(time_idx_ofinterest(maxloc))]; %maxv',%time_points(time_idx_ofinterest(maxloc))'
        end
        for ic = odor_picks{iodor}
            idxC = find(synch_condition{mouse_i}(:,1) == concentrations(ic)); 
            idx = intersect(idxO, idxC); 

            intensity_class = [intensity_class; ic-3]; 
            
            avg_profile = mean(synchrony_timeseries{mouse_i}(idx, :), 1); 
            [maxv, maxloc] = max(avg_profile(time_idx_ofinterest)');
            synch_class{1} = [synch_class{1};  maxv, time_points(time_idx_ofinterest(maxloc))]; %maxv',%time_points(time_idx_ofinterest(maxloc))'
        end
    end
end
ppm_concentration_class = zeros(size(ppm_concentration_values)); 
ppm_concentration_class(ppm_concentration_values<=t1) = 1; 
ppm_concentration_class(ppm_concentration_values>t2) = 3; 
ppm_concentration_class(intersect(find(ppm_concentration_values>t1), find(ppm_concentration_values<=t2))) = 2; 

figure
subplot(1, 2, 1)
hist(intensity_class)
ylim([0, 50])
subplot(1, 2, 2)
hist(ppm_concentration_class)
ylim([0, 50])


%% Apply neural net classifier

Y = {intensity_class, ppm_concentration_class};     % 1 x N
%Y_target = ppm_concentration_class ; 
Nneurons = 30; 

for ioutput = 1 : size(Y, 2)
    % PERCEPTRON
    X_input = synch_class{ioutput}';      % 2 x N
    Y_target = Y{ioutput};     % 1 x N
    %Y_target = ppm_concentration_class ; 
    T = full(ind2vec(Y_target'));
    % SVM
    %X_input = synch_class{ioutput};      % 2 x N
    %Y_target = Y{ioutput};     % 1 x N
    %Y_target = ppm_concentration_class ; 
    %T = full(ind2vec(Y_target'))';
    for irun = 1 : 10
       
        % PERCEPTRON
        clear net
        % Create a feedforward network with 10 hidden neurons
        net = patternnet(Nneurons);  % 10 hidden neurons (you can tune this)
        % Optional: Set transfer functions
        net.layers{1}.transferFcn = 'tansig';    % hidden layer
        net.layers{2}.transferFcn = 'softmax';   % output layer (for probabilities
        
        % Optional: Divide data
        net.divideParam.trainRatio = 0.8;
        net.divideParam.valRatio   = 0.1;
        net.divideParam.testRatio  = 0.1;
        
        % Train the network
        [net, tr] = train(net, X_input, T);  % X' because patternnet expects [features × samples]
        
        % Predict on training data
        Y_pred = net(X_input);
        [~, idx] = max(Y_pred,[],  1);          % Find index of max
        Y_out = zeros(size(Y_pred));  
        for i = 1:length(idx)% Initialize with zeros
            Y_out(idx(i), i) = 1;    
        end
        
        [~, class_pred] = max(Y_pred);   % Get class indices
        
        % Optional: Confusion matrix
        %plotconfusion(T, Y_pred);
        
        %Optional: Classification accuracy
        true_labels = vec2ind(T);  % Convert one-hot back to class labels
        pred_labels = vec2ind(Y_out); 
        accuracy(ioutput, irun) =sum(class_pred == true_labels) / length(true_labels);
        fprintf('Accuracy: %.2f%%\n', accuracy * 100);
        CM = confusionmat(pred_labels,true_labels); 
        CM_percentage = CM./sum(CM, 1);
        
        allCM{ioutput}(:,:,irun) = CM_percentage; 
    end
end
%
figure
for ioutput = 1 : length(allCM)
    subplot(1, length(allCM), ioutput)
    mycm = mean(allCM{ioutput}, 3); 
    CM_percentage = mycm./sum(mycm, 1); 
    imagesc(CM_percentage, [0, 1])
    hold on
    for i = 1:size(CM_percentage, 1)
        for j = 1:size(CM_percentage, 2)
            value = CM_percentage(i,j);
            text(j, i, num2str(value), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontWeight', 'bold', ...
                'Color', 'black');  % You can change to 'black' if needed
        end
    end
end
colormap("copper")
figure, boxplot(accuracy')
title(num2str(ranksum(accuracy(1,:), accuracy(2,:))))

%% Plot dependency of accuracy from Nneurons
ineuron = 0; 
temp_result = []; 
nneurons_accuracy = []
xneurons = 4:2:50; 
ioutput = 1; 
for Nneurons = xneurons
    ineuron = ineuron +1; 
    X_input = synch_class{ioutput}';      % 2 x N
    Y_target = Y{ioutput};     % 1 x N
    %Y_target = ppm_concentration_class ; 
    T = full(ind2vec(Y_target'));
    for irun = 1 : 5
    % PERCEPTRON
        clear net
        % Create a feedforward network with 10 hidden neurons
        net = patternnet(Nneurons);  % 10 hidden neurons (you can tune this)
        % Optional: Set transfer functions
        net.layers{1}.transferFcn = 'tansig';    % hidden layer
        net.layers{2}.transferFcn = 'softmax';   % output layer (for probabilities
        
        % Optional: Divide data
        net.divideParam.trainRatio = 0.8;
        net.divideParam.valRatio   = 0.1;
        net.divideParam.testRatio  = 0.1;
        
        % Train the network
        [net, tr] = train(net, X_input, T);  % X' because patternnet expects [features × samples]
        
        % Predict on training data
        Y_pred = net(X_input);
        [~, idx] = max(Y_pred,[],  1);          % Find index of max
        Y_out = zeros(size(Y_pred));  
        for i = 1:length(idx)% Initialize with zeros
            Y_out(idx(i), i) = 1;    
        end
        
        [~, class_pred] = max(Y_pred);   % Get class indices
        
        % Optional: Confusion matrix
        %plotconfusion(T, Y_pred);
        
        %Optional: Classification accuracy
        true_labels = vec2ind(T);  % Convert one-hot back to class labels
        pred_labels = vec2ind(Y_out); 
        temp_result(irun) =sum(class_pred == true_labels) / length(true_labels);
    end
    nneurons_accuracy(ineuron)= mean(temp_result); 
end

figure
plot(xneurons, nneurons_accuracy, 'o-')


