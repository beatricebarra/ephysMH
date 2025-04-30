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
for ifilename = 1 : length(filenames)
     
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
    good_units_PCx{ifilename} = PC_units_idx; 
end

%% Read complete dataset with all mice
savefilenameallmice = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/all_mice_dataset_OB.mat'; 
load(savefilenameallmice)
sig = 0.010; 
mint = 0; 
maxt = 0.5;
tvect = linspace(mint, maxt, (maxt-mint)/sig); 

%% Try out your new synch measure on each mouse and save results

mycolormap = jet(10); 
mintime = -0.2; 
maxtime = 0.5; 
results = {};
synchrony_timeseries = {}; 
synch_condition = {}; 

dt = 0.001; 
maxW = 0.01; 
synch_window =0.2; 
for mouse_i= 1:7
    iunits = find(imouse == mouse_i); 
    % Find which units are responsive for each odor
    [selunits] = select_responsive_units(all_mice_spikes_table, all_mice_condition_table, iunits); 
    % index of conditions of this mouse (they are all the same technically
    iconditions = (mouse_i-1)*32+1:(mouse_i)*32;
    % Initializations 
    results{mouse_i} = []; 
    synchrony_timeseries{mouse_i} = []; 
    synch_condition{mouse_i} = [];  
    
    figure
    irow = 0;
    irow_rep = 0; 
    for itrial = 1 : size(all_mice_spikes_table, 2) % this is the condition, i.e all trial types % 1 to 32
        % Finding the minimum number of repetitions for this trial (some units have less somehow, I think this is a bug
        Nrep = length(all_mice_spikes_table{1, itrial}); 
        for iu = 1 : size(all_mice_spikes_table, 1)
            Nrep = min([Nrep, length(all_mice_spikes_table{iu, itrial})]); 
        end
        % Take odor index of this trial type
        odor = all_mice_condition_table(iconditions(itrial), 1); 
        % For each repetition
        sychvect = []; 
        for irep = 1 : Nrep
            disp(['Trial ', num2str(itrial), ' Rep = ',  num2str(irep)])
            
            % Selecting units of interest. Here I am merging all odors right
            % now, but I could filter for units responsive to this odor
            thistrialtable = {};
            ii= 0; 
            for iu = selunits{odor}
                ii = ii+1; 
                idx = intersect(find(all_mice_spikes_table{iu, itrial}{irep}>mintime), find(all_mice_spikes_table{iu, itrial}{irep}<maxtime)); 
                thistrialtable{ii} = all_mice_spikes_table{iu, itrial}{irep}(idx); 
                
            end
            
            time_points = [mintime:dt:maxtime]; 
            total = floor((maxtime-mintime)/dt); 
            checkpoints = floor(linspace(0.01, 1, 100) * total);  % 10%, 20%, ..., 100%    
            % Compute synchrony
            [Synch] = computeSynchMeasure(thistrialtable, time_points, maxW); 
            % Take max of synch
            synchpeak = max(Synch(intersect(find(time_points > 0), find(time_points < synch_window)))); 
            % Save time profiles of synchrony measure
            irow_rep = irow_rep +1; 
            synchrony_timeseries{mouse_i}(irow_rep,:) = Synch; 
            synch_condition{mouse_i}(irow_rep,:) =  [all_mice_condition_table(iconditions(itrial), 3), all_mice_condition_table(iconditions(itrial), 1)]; 

            % save features for future plotting -- 
            %if I do it like this,each repetition is a data point
            %irow = irow+1; 
            %results{mouse_i}(irow,:) = [synchpeak, all_mice_condition_table(iconditions(itrial), 3), all_mice_condition_table(iconditions(itrial), 1)]; 
            % Otherwise .. 
            sychvect = [sychvect; synchpeak]; 

        end
        %if I do it like this,each repetition is a data point
        irow = irow+1; 
        results{mouse_i}(irow,:) = [mean(sychvect), all_mice_condition_table(iconditions(itrial), 3), all_mice_condition_table(iconditions(itrial), 1)]; 
        
    

    end
    
    % Plot changes in synchrony with concentration
    figure
    hold on
    mycolormap = lines(5); 
    theodors = unique(results{mouse_i}(:,3)); 
    for iodor = 1 : length(unique(theodors))
        idxO = find(results{mouse_i}(:,3) == theodors(iodor)); 
        concentrations = unique(results{mouse_i}(idxO,2)); 
        means = []; 
        for ic = 1 : length(concentrations)
            idxC = find(results{mouse_i}(:,2) == concentrations(ic)); 
            idx = intersect(idxO, idxC); 
            %plot((results(idx,2)), (results(idx, 1)),  'o', 'Color', mycolormap(iodor,:))
            means = [means; mean(results{mouse_i}(idx, 1))]; 
    
            errorbar(mean(results{mouse_i}(idx,2)), mean(results{mouse_i}(idx, 1)), std(results{mouse_i}(idx, 1))/sqrt(length(results{mouse_i}(idx, 1))), ...
                '-o', 'Color', mycolormap(iodor,:), 'MarkerFaceColor', mycolormap(iodor,:))
        end
        plot(concentrations, means, '-o', 'Color', mycolormap(iodor,:), 'MarkerFaceColor', mycolormap(iodor,:))
    end
    set(gca, 'XScale', 'log')
    xlabel('Concentration')
    ylabel('Synchrony')

end

%% Plot synchrony time series for different mice 
mycolormaps = generateColorMaps(10); 
colorlabel = {'red', 'blue', 'yellow', 'green'}; 
figure
for mouse_i = 1: 7
    for iodor = 1 : length(unique(theodors))
        idxO = find(synch_condition{mouse_i}(:,2) == theodors(iodor)); 
        concentrations = unique(synch_condition{mouse_i}(idxO,1)); 
        for ic = 1 : length(concentrations)
            idxC = find(synch_condition{mouse_i}(:,1) == concentrations(ic)); 
            idx = intersect(idxO, idxC); 
            ax(mouse_i, iodor) = subplot(7, length(unique(theodors)), (mouse_i-1)*4 + iodor); 
            hold on
            if length(idx)> 1
                stdshade(synchrony_timeseries{mouse_i}(idx, :), 0.3, mycolormaps.(colorlabel{iodor})(ic,:))
            end
            avg_profile = mean(synchrony_timeseries{mouse_i}(idx, :))
        end
    end
end
linkaxes(ax, 'x')

%% Create dataset for regressor that includes peak value and location of synchrony time_series
% Each repetition counts as a separate time point
odor_symbols = ['o', 's', 't', 'x']; 
intensity_file  = '/Users/barrab01/Documents/Repos/ephysMH/perceived_intensity_curves.xlsx'; 
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
for mouse_i = 1 : 7
    odor_concentrations = fliplr(synch_condition{mouse_i}); 
    realI = computeIntensity(intensity_file, odor_concentrations); 
    
    
    mouse_label = [mouse_label; ones(size(synch_condition{mouse_i}(:,1)))]; 
    odor_label = [odor_label; odor_concentrations(:,1)]; 
    concentration = [concentration; odor_concentrations(:, 2)]; 
    temp = []; 
    for ii = 1:length(odor_concentrations(:,1))
        io = odor_concentrations(ii,1); 
        temp(ii) = odor_concentrations(ii, 2)*ppm_col(io); 
    end
    ppm_concentration = [ppm_concentration; temp]; 
    intensity = [intensity; realI']; 
    [maxv, maxloc] = max(synchrony_timeseries{mouse_i}(:, time_idx_ofinterest)')
    synch = [synch; maxv', time_points(time_idx_ofinterest(maxloc))']; 
end

       
idx_nonnegative = find(concentration~=0); 
Y = log10(concentration(idx_nonnegative)); 
Y = (ppm_concentration(idx_nonnegative)); 
Y = intensity(idx_nonnegative); 
X = [synch(idx_nonnegative, :);]; 
%X = [gfeats(idx_nonnegative,:)];

figure
for i = 1 : size(X, 2)
    subplot(1, 2, i)
    hold on
    plot(X(:,i), Y, 'o')
    linearCoef = polyfit(X(:, i),Y,2);
    linearFit = polyval(linearCoef,X(:, i));
    hold on
    plot(X(:, i),linearFit,'k-', 'Linewidth', 2)
end


figure;
hold on; 
MSE_realInt = []; 
for i = 1 : 5
    [rfModel, X_test, Y_test, Y_pred, mse] = emh_random_forest_regressor(X, Y, odor_label, odor_symbols, '-'); 
    MSE_realInt(i) = mse; 
end

%% Create dataset for regressor that includes peak value and location of synchrony time_series
% Repetitions are averaged to create a single datapoint
odor_symbols = ['o', 's', 't', 'x']; 
intensity_file  = '/Users/barrab01/Documents/Repos/ephysMH/perceived_intensity_curves.xlsx'; 
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
for mouse_i = 1: 7
    for iodor = 1 : length(unique(theodors))
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


%% Create intensity regression model for the 2 features separately 
figure
idx_nonnegative = intersect(find(concentration~=0), find(intensity~=0)); 
for i = 1 : size(synch, 2)
    subplot(1, 2, i)
    hold on
    x = synch(idx_nonnegative,i); 
    y = log10(intensity(idx_nonnegative)); 
    plot(x, y, 'ro')
    %plot(synch(:,i), ppm_concentration, 'bo')
    linearCoef_intensitymodel(i, :) = polyfit(x,y,1);
    linearFit = polyval(linearCoef_intensitymodel(i,:),x);
    hold on
    plot(x,linearFit,'k-', 'Linewidth', 2)
    title(['y = ', num2str(linearCoef_intensitymodel(i,1)), '*x + ', num2str(linearCoef_intensitymodel(i, 2))])
end
I = log10(intensity(idx_nonnegative)); 
Icap = zeros(size(idx_nonnegative, 1), 1); 
for i = 1 : size(synch, 2)
    Icap = Icap + linearCoef_intensitymodel(i, 1).*synch(idx_nonnegative, i) + linearCoef_intensitymodel(i, 2) 
end
Icap = Icap./2; 
figure
plot(I, Icap, 'o')
hold on
plot([-10, 1], [-10, 1], 'r--')
lC = polyfit(I,Icap,1);
linearFit = polyval(lC,I);
plot(I, linearFit, 'k-')
rmse = sqrt(median((I - Icap).^2));
nrmse = rmse/(max([I]) - min([Icap]));

%% If I want to average measures 
clear net
close all

idx_nonnegative = (find(concentration~=0)); 
Y = [log10(intensity(idx_nonnegative)), (log10(concentration(idx_nonnegative))),  (log10(ppm_concentration(idx_nonnegative)))]; 
X = [synch(idx_nonnegative, :)]; 

for ioutput = 1 : size(Y, 2)

    figure;
    % PERCEPTRON
    %X_input = X';      % 2 x N
    %Y_target = Y(:, ioutput)';     % 1 x N
    
    % SVM
    X_input = X;      % 2 x N
    Y_target = Y(:, ioutput);     % 1 x N
    
    %RANDOM FOREST
    %X_input = X;      % 2 x N
    %Y_tar = Y(:, ioutput); 

    for i =  1 : 20
        % PERCEPTRON
%         Nneurons  = floor(size(X, 1)/5); 
%         Nneurons = 100; 
%         % Create a feedforward network with 10 hidden neurons
%         net = fitnet(Nneurons);  % You can try other sizes too
%         net.trainParam.max_fail = 10;
%         net.layers{1}.transferFcn = 'tansig';  % (default) Hidden layer
%         net.layers{2}.transferFcn = 'purelin';%purelin
%         
%         % Optional: split data manually (default is 70/15/15 train/val/test)
%         net.divideParam.trainRatio = 0.8;
%         net.divideParam.valRatio = 0.1;
%         net.divideParam.testRatio = 0.1;
%         
%         % Train the network
%         [net, tr] = train(net, X_input, Y_target);
%         
%         % Predict on training data
%         Y_pred = net(X_input);
        

        % SVM
        %svmModel = fitrsvm(X_input, Y_target);  % uses default settings (linear kernel, etc.)
        trainRatio = 0.8;
        idx = randperm(length(Y_target));
        trainIdx = idx(1:round(trainRatio * length(Y_target)));
        testIdx = idx(round(trainRatio * length(Y_target)) + 1:end);
        
        X_train = X_input(trainIdx, :);
        Y_train = Y_target(trainIdx);
        X_test = X_input(testIdx, :);
        Y_test = Y_target(testIdx);
        svmModel = fitrsvm(X_input, Y_target, ...
            'KernelFunction', 'polynomial', ...
            'KernelScale', 'auto', ...
            'Standardize', true);
        Y_pred = predict(svmModel, X_input);
        

        % RANDOM FOREST
        %[rfModel, X_test, Y_target, Y_pred, nrmse] = emh_random_forest_regressor(X_input, Y_tar, odor_label, odor_symbols, '-'); 

        % Plot actual vs predicted
        plot(Y_target, Y_pred, 'o');
        xlabel('True Y'); ylabel('Predicted Y');
        title('Regression: True vs Predicted');
        grid on; axis equal;
        hold on
        plot(Y_target, Y_target, 'r', 'Linewidth', 2)
        
        linearCoef = polyfit(Y_target,Y_pred,1);
        linearFit = polyval(linearCoef,Y_target);
        hold on
        plot(Y_target,linearFit,'k-', 'Linewidth', 2)
        
        % Compute mean squared error
        mse_error = mse(Y_pred - Y_target);
        
        rmse = sqrt(median((Y_target - Y_pred).^2));
        nrmse = rmse/(max([Y_target]) - min([Y_target]));
        %mape = mean(abs((Y_test - Y_pred) ./ Y_test)) * 100;
        net_errors(i , ioutput) = nrmse; 
        net_rmse(i,ioutput ) = rmse; 
        
        R = corrcoef(Y_pred, Y_target);
        Rcoeff(i , ioutput)  = R(1,2);  % This is the Pearson correlation coefficient
        %fprintf('Mean Squared Error: %.4f\n', nrmse);
    end
end
%fprintf('Mean Squared Error: %.4f\n', mse_error);
figure
subplot(1, 2, 1)
boxplot(net_errors)
subplot(1, 2, 2)
boxplot(Rcoeff)


%% Create dataset for CLASSIFIER that includes peak value and location of synchrony time_series
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
        
        for ic = 4 : 6 % defin concentration classes 
            temp = [temp; log10(concentrations(ic)*ppm_col(iodor))]; 
        end
    end
   
    
end  
step = (1.1*max(temp) - 1.1*min(temp)) / 3;
t1 = min(temp) + step;      % Lower threshold
t2 = min(temp) + 2 * step; 
figure, hold on, plot(temp, rand(size(temp)), 'o'), plot([t1, t1], [0, 1]), plot([t2, t2], [0, 1])
plot([min(temp), min(temp)], [0, 1])

odor_picks_ppm = {[3 : 6], [3 : 6], [3 : 5], [3 : 6]}; 
odor_picks = {[4 : 6], [4 : 6], [4 : 5], [4 : 6]}; 

for mouse_i = 1: 7    % Upper threshold
    for iodor = 1 : length(unique(theodors))
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


Y = {intensity_class, ppm_concentration_class};     % 1 x N
%Y_target = ppm_concentration_class ; 
Nneurons = 5; 

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
    for irun = 1 : 20
       
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
        plotconfusion(T, Y_pred);
        %Optional: Classification accuracy
        true_labels = vec2ind(T);  % Convert one-hot back to class labels
        pred_labels = vec2ind(Y_out); 
        accuracy(ioutput, irun) =sum(class_pred == true_labels) / length(true_labels);
        fprintf('Accuracy: %.2f%%\n', accuracy * 100);
        CM = confusionmat(pred_labels,true_labels); 
        CM_percentage = CM./sum(CM, 1); 


        % SVM 

%         trainRatio = 0.8;
%         idx = randperm(length(Y_target));
%         trainIdx = idx(1:round(trainRatio * length(Y_target)));
%         testIdx = idx(round(trainRatio * length(Y_target)) + 1:end);
%         
%         X_train = X_input(trainIdx, :);
%         Y_train = Y_target(trainIdx);
%         X_test = X_input(testIdx, :);
%         Y_test = Y_target(testIdx);
%         svmModel = fitcecoc(X_train, Y_train);
% 
%         Y_pred = predict(svmModel, X_test);
% 
%        
%         
%         acc_svm = mean(Y_pred == Y_test);
%         accuracy(ioutput, irun) = mean(Y_pred == Y_test);
%         fprintf('Accuracy: %.2f%%\n', accuracy * 100);
%         CM = confusionmat(Y_pred,Y_test); 
%         CM_percentage = CM./sum(CM, 1); 
        
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
figure, boxplot(accuracy')
title(num2str(ranksum(accuracy(1,:), accuracy(2,:))))
%% Plot mouse specific results
fig1 = figure; 
fig2 = figure; 
fig3 = figure; 
for mouse_i = 1: 7
    % Plot changes in synchrony with concentration
    figure
    hold on
    mycolormap = lines(5); 
    theodors = unique(results{mouse_i}(:,3)); 
    for iodor = 1 : length(unique(theodors))
        idxO = find(results{mouse_i}(:,3) == theodors(iodor)); 
        concentrations = unique(results{mouse_i}(idxO,2)); 
        means = []; 
        for ic = 1 : length(concentrations)
            idxC = find(results{mouse_i}(:,2) == concentrations(ic)); 
            idx = intersect(idxO, idxC); 
            %plot((results(idx,2)), (results(idx, 1)),  'o', 'Color', mycolormap(iodor,:))
            means = [means; mean(results{mouse_i}(idx, 1))]; 
    
            errorbar(mean(results{mouse_i}(idx,2)), mean(results{mouse_i}(idx, 1)), std(results{mouse_i}(idx, 1))/sqrt(length(results{mouse_i}(idx, 1))), ...
                '-o', 'Color', mycolormap(iodor,:), 'MarkerFaceColor', mycolormap(iodor,:))
        end
        plot(concentrations, means, '-o', 'Color', mycolormap(iodor,:), 'MarkerFaceColor', mycolormap(iodor,:))
    end
    set(gca, 'XScale', 'log')
    xlabel('Concentration')
    ylabel('Synchrony')
    figure(fig1)
    for iodor = 1 : length(unique(theodors))
        subplot(1, 4, iodor)
        hold on
        idxO = find(results{mouse_i}(:,3) == theodors(iodor)); 
        concentrations = unique(results{mouse_i}(idxO,2)); 
        means = []; 
        for ic = 1 : length(concentrations)
            idxC = find(results{mouse_i}(:,2) == concentrations(ic)); 
            idx = intersect(idxO, idxC); 
            %plot((results(idx,2)), (results(idx, 1)),  'o', 'Color', mycolormap(iodor,:))
            means = [means; mean(results{mouse_i}(idx, 1))]; 
    
            errorbar(mean(results{mouse_i}(idx,2)), mean(results{mouse_i}(idx, 1)), std(results{mouse_i}(idx, 1))/sqrt(length(results{mouse_i}(idx, 1))), ...
                '-o', 'Color', mycolormap(iodor,:), 'MarkerFaceColor', mycolormap(iodor,:))
        end
        plot(concentrations, means, '-o', 'Color', mycolormap(iodor,:), 'MarkerFaceColor', mycolormap(iodor,:))
        set(gca, 'XScale', 'log')
        xlabel('Concentration')
        ylabel('Synchrony')
    end

    figure(fig2)
    for iodor = 1 : length(unique(theodors))
        subplot(1, 4, iodor)
        hold on
        idxO = find(results{mouse_i}(:,3) == theodors(iodor)); 
        x = log10(results{mouse_i}(idxO, 2)); 
        y = results{mouse_i}(idxO, 1); 
        idx = ~isinf(x); 
        x = x(idx); 
        y = y(idx);
        plot(x, y, 'o')
        linearCoef = polyfit(x,y,1);
        linearFit = polyval(linearCoef,x);
        plot(x,y,'s', x,linearFit,'r-')
        %set(gca, 'XScale', 'log')
        xlabel('Concentration')
        ylabel('Synchrony')
    end
    
end


allmice_means = []; 
for iodor = 1 : length(unique(theodors))
    allX =  [];
    allY = [];
    
    hold on
    for mouse_i = 1: 7
        
        idxO = find(results{mouse_i}(:,3) == theodors(iodor)); 
        x = log10(results{mouse_i}(idxO, 2)); 
        y = results{mouse_i}(idxO, 1); 
        idx = ~isinf(x); 
        x = x(idx); 
        y = y(idx);
        allX = [allX, x]; 
        allY = [allY, y]; 

        concentrations = unique(results{mouse_i}(idxO,2)); 
        means = []; 
        for ic = 1 : length(concentrations)
            idxC = find(results{mouse_i}(:,2) == concentrations(ic)); 
            idxOC = intersect(idxO, idxC); 
            %plot((results(idx,2)), (results(idx, 1)),  'o', 'Color', mycolormap(iodor,:))
            means = [means; mean(results{mouse_i}(idxOC, 1))]; 
        end
        allmice_means = [allmice_means; means']; 
    end
    % Figure with all datapoints and regression
    figure(fig2)
    subplot(1, 4, iodor)
    linearCoef = polyfit(allX,allY,1);
    linearFit = polyval(linearCoef,allX);
    plot(allX,allY,'s', allX,linearFit,'k-', 'Linewidth', 2)
    xlabel('Concentration')
    ylabel('Synchrony')
    figure(fig3)
    errorbar(log10(concentrations), mean(allmice_means), std(allmice_means), ...
                '-o', 'Color', mycolormap(iodor,:), 'MarkerFaceColor', mycolormap(iodor,:), 'Markersize', 10, 'Linewidth', 2)
    %set(gca, 'XScale', 'log')
    xlabel('Concentration')
    ylabel('Synchrony')
end
%% Plot boxplots of synch at intensity -matched concentrations
% Initialize
boxplot_data = {}; 
for i = 1:3
    boxplot_data{i} = []; 
    group_label{i} = []; 
end

for iodor = 1 : length(unique(theodors))
    for mouse_i = 1: 7
        idxO = find(results{mouse_i}(:,3) == theodors(iodor)); 
        concentrations = unique(results{mouse_i}(idxO,2)); 
        iconc = 1; 
        for ic = 4:6
            idxC = find(results{mouse_i}(:,2) == concentrations(ic)); 
            idxOC = intersect(idxO, idxC); 
            boxplot_data{iconc } = [boxplot_data{iconc}; results{mouse_i}(idxOC,1)]; 
            group_label{iconc } = [group_label{iconc}; iodor*ones(size(results{mouse_i}(idxOC,1)))]; 
            iconc = iconc+1; 

        end
    end
end

figure
hold on
for ic = 1:3
    subplot(1, 3, ic)
    hold on
    iodor
    boxplot(boxplot_data{ic},group_label{ic})
    ylim([0, 0.4])
end
% Linearize data for statistical comparison
stats_data = []; 
stats_group = []; 

for ic = 1:3
   stats_data = [stats_data; boxplot_data{ic}]; 
   index = (ic-1)*length(unique(group_label{ic})) + group_label{ic}; 
   stats_group = [stats_group; index]; 
end

[P,ANOVATAB,STATS] = kruskalwallis(stats_data,stats_group); 
multcompare(STATS)

%% Plot results as synch = f(intensity)

% Add column to "results" with intensity label
intensity_file  = '/Users/barrab01/Documents/Repos/ephysMH/perceived_intensity_curves.xlsx'; 
ppm_col = [5617.105263; 16842.10 ; 526.16; 5065.789474];
for mouse_i = 1 : length(results)
    odor_concentrations = fliplr(results{mouse_i}(:, 2:3)); 
    realI = computeIntensity(intensity_file, odor_concentrations); 
    results{mouse_i}(:, 4) = realI; 

    plot(results{mouse_i}(:, 4), results{mouse_i}(:, 1), 'o')
end
figure
hold on
for mouse_i = 1 : length(results)
    plot((results{mouse_i}(:, 4)), results{mouse_i}(:, 1), 'o')

end



%% Train random forest regressor on classifying concentration and intensity 
odor_symbols = ['o', 's', 't', 'x']; 
intensity_file  = '/Users/barrab01/Documents/Repos/ephysMH/perceived_intensity_curves.xlsx'; 
ppm_col = [5617.105263; 16842.10 ; 526.16; 5065.789474];
X = []; 
Y = []; 
mouse_label = []; 
odor_label = []; 
intensity = []; 
concentration = []; 
synch = []; 
ppm_concentration = []; 
for mouse_i = 1 : 7
    odor_concentrations = fliplr(results{mouse_i}(:, 2:3)); 
    realI = computeIntensity(intensity_file, odor_concentrations); 
    
    
    mouse_label = [mouse_label; ones(size(results{mouse_i}(:,1)))]; 
    odor_label = [odor_label; odor_concentrations(:,1)]; 
    concentration = [concentration; odor_concentrations(:, 2)]; 
    temp = []; 
    for ii = 1:length(odor_concentrations(:,1))
        io = odor_concentrations(ii,1); 
        temp(ii) = odor_concentrations(ii, 2)*ppm_col(io); 
    end
    ppm_concentration = [ppm_concentration; temp]; 
    intensity = [intensity; realI']; 
    synch = [synch;results{mouse_i}(:,1) ]; 
end

       
idx_nonnegative = find(concentration~=0); 
Y = concentration(idx_nonnegative); 
Y = ppm_concentration(idx_nonnegative); 
Y = intensity(idx_nonnegative); 
X = [synch(idx_nonnegative);]; 
%X = [gfeats(idx_nonnegative,:)];


figure;
hold on; 
MSE_realInt = []; 
for i = 1 : 5
    [rfModel, X_test, Y_test, Y_pred, mse] = emh_random_forest_regressor(X, Y, odor_label, odor_symbols, '-'); 
    MSE_realInt(i) = mse; 
end


%% END OF PROGRESS
% old from now on
%% NOW DO THE Same BUT MIX UNITS FROM ALL MICE ALL WHILE FILTERING THE UNITS THAT ANSWER TO SPECIFIC ODOR

dt = 0.001; %median(diff(sort(multivariate_spiketrain))); 
maxW = 0.01; 

irow = 0; 
% Find units that are responsive across at least one condition of an
% odorant
responsiveness = []
condition_table = all_mice_condition_table(1:32,:); 
odors = unique(condition_table(:,1)); 
selunits = {}; 
for iodor = 1 : length(unique(condition_table(:,1)))
    selunits{iodor} = []; 
end
for iu = 1 : size(all_mice_spikes_table, 1)
    disp(['Unit ', num2str(iu)])
    for iodor = 1 : length(unique(condition_table(:,1)))
        idxO = find(condition_table(:,1) == odors(iodor)); 
        concs = unique(condition_table(idxO,3));
        concs = concs(concs>0); 
        for iconc = 1: length(concs)
            idxC = find(condition_table(:,3) == concs(iconc)); 
            idxT = intersect(idxO, idxC); 
            
            % Compute whether unit is responsive in this condition
            [response_index, pvalue] = find_responsive_units( all_mice_spikes_table{iu, idxT}, [0, 0.15], [-0.2, -0.05]); 
            responsiveness(iu, iodor, iconc) = pvalue;

        end
        % Find if there is a response for any of the concentrations on
        % this odorant
        if any(responsiveness(iu, iodor, :)<0.05) % save unit as viable only if responsive
            selunits{iodor} = [selunits{iodor}, iu] 
        end
    end
end
%% Plot repartition of all responsive units, by odor
all_units = [1 : size(all_mice_spikes_table, 1)]
sometimes_responsive = selunits{1}; 
for iodor = 2 : 4
    sometimes_responsive =union(sometimes_responsive, selunits{iodor});  
end

never_responsive = length(all_units) - length(sometimes_responsive); 
multiple_responses = length(all_units) - never_responsive;
for iodor = 1 : 4
    exclusive_responses{iodor} = selunits{iodor}; 
    other_odors =setdiff([1:1:4], [iodor]); 
    for oodor = other_odors
        exclusive_responses{iodor} =setdiff(exclusive_responses{iodor}, selunits{oodor});  
    end
    multiple_responses = multiple_responses-length(exclusive_responses{iodor}); 
end


figure
pie([ never_responsive, multiple_responses,...
    length(exclusive_responses{1}), length(exclusive_responses{2}), ...
    length(exclusive_responses{3}), length(exclusive_responses{4})])
pielabels = {'No resp', 'Many', 'odor1', 'odor2', 'odor3', 'odor4'}
legend(pielabels)
sum([multiple_responses, never_responsive, length(exclusive_responses{1}), length(exclusive_responses{2}), length(exclusive_responses{3}), length(exclusive_responses{4})])

%% Compute synchrony for cells specific odorants, across all mice
irow = 0; 
mintime = -0.2; 
maxtime = 0.5; 
for itrial = 1 : size(all_mice_spikes_table, 2)
    
    Nrep = length(all_mice_spikes_table{1, itrial}); % Finding the minimum number of repetitions for this trial (some units have less somehow, I think this is a bug
    for iu = 1 : size(all_mice_spikes_table, 1)
        Nrep = min([Nrep, length(all_mice_spikes_table{iu, itrial})]); 
    end
    iodor = condition_table(itrial, 1); 
    % For each repetition
    for irep = 1 : Nrep
        disp(['Trial ', num2str(itrial), ' Rep = ',  num2str(irep)])
        irow = irow+1; 
        ii=0; 
        thistrialtable = {};
        for iu = selunits{iodor}
            ii = ii+1; 
            idx = intersect(find(all_mice_spikes_table{iu, itrial}{irep}>mintime), find(all_mice_spikes_table{iu, itrial}{irep}<maxtime)); 
            thistrialtable{ii} = all_mice_spikes_table{iu, itrial}{irep}(idx); 
        end
        
        time_points = [mintime:dt:maxtime]; 
        % Compute synchrony
        [Synch] = computeSynchMeasure(thistrialtable, time_points, maxW); 
        % Take max of synch
        synchpeak = max(Synch(intersect(find(time_points > 0), find(time_points < 0.2)))); 
        % save features for future plotting
        results_across_mice(irow,:) = [synchpeak, condition_table(itrial, 3), condition_table(itrial, 1)]; 
        
    end
end
%% Plot results of synchrony computed on all mice together
figure
theodors = unique(results_across_mice(:,3)); 
for iodor = 1 : length(unique(results_across_mice(:,3)))
    allX =  [];
    allY = [];

    idxO = find(results_across_mice(:,3) == theodors(iodor)); 
    x = log10(results_across_mice(idxO, 2)); 
    y = results_across_mice(idxO, 1); 
    idx = ~isinf(x); 
    allX = x(idx); 
    allY = y(idx);
    
    concentrations = unique(results_across_mice(idxO,2)); 
    means = []; 
    for ic = 1 : length(concentrations)
        idxC = find(results_across_mice(:,2) == concentrations(ic)); 
        idxOC = intersect(idxO, idxC); 
        %plot((results(idx,2)), (results(idx, 1)),  'o', 'Color', mycolormap(iodor,:))
        means = [means; mean(results_across_mice(idxOC, 1))]; 
    end
    
    subplot(1, 4, iodor)
    linearCoef = polyfit(allX,allY,1);
    linearFit = polyval(linearCoef,allX);
    hold on
    plot(allX,allY,'s', allX,linearFit,'k-', 'Linewidth', 2)
    plot(log10(concentrations),means,'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 10 )
    xlabel('Concentration')
    ylabel('Synchrony')
   
 
end



%% ALSO DO PCA

%%

[value_max, idx_closest_spike] = cellfun(@(x) min(abs(x-t)) , thistrialtable , 'UniformOutput' , false); 
    
% Pre-allocate new_array once
new_array = cell(size(idx_closest_spike));

% Logical index of empty cells
empty_idx = cellfun(@isempty, idx_closest_spike);

% Fill empty cells with maxW
new_array(empty_idx) = {maxW};

% Fill non-empty cells with min(value_max, maxW)
nonempty_idx = ~empty_idx;
clipped_values = cellfun(@(v) min(v, maxW), value_max(nonempty_idx), 'UniformOutput', false);
new_array(nonempty_idx) = clipped_values;

%% Take a trial and test the spike measure
multivariate_spiketrain = []
thistrialtable = {}

mintime = -1; 
maxtime = 1; 

dt = 0.001; %median(diff(sort(multivariate_spiketrain))); 

itrial = 6; 
for iu = 1 : size(all_mice_spikes_table, 1)
    multivariate_spiketrain = [multivariate_spiketrain; all_mice_spikes_table{iu, itrial}{1}]; 
    idx = intersect(find(all_mice_spikes_table{iu, itrial}{1}>mintime), find(all_mice_spikes_table{iu, itrial}{1}<maxtime)); 
    thistrialtable{iu} = all_mice_spikes_table{iu, itrial}{1}(idx); 
end
figure
hist(diff(sort(multivariate_spiketrain)), 1000)


total = floor((maxtime-mintime)/dt); 
checkpoints = floor(linspace(0.01, 1, 100) * total);  % 10%, 20%, ..., 100%
SPIKE = []
tvect = []
time_points = [mintime:dt:maxtime]; 
%[D_Sm, S_m_t] = multivariate_average_spike_distance(thistrialtable, max(timepoints) - min(time_points), time_points)

%S = multivariateSPIKEDistance(thistrialtable, time_points)
far_points = [mintime*10, maxtime*10]; 
for it = 1 : length(time_points) % for each time instant

    t = time_points(it ); % Compute time 
    if any(it == checkpoints)
        fprintf('Progress: %.0f%% complete\n', (it / total) * 100);
        pause(0.01)
    end
    tvect(it) = t; 
    %cellfun(@(x) find(x<t) , thistrialtable)
    %Find preceding spikes 
    preceding_spikes = cellfun(@(x) x(find(x<t)) , thistrialtable , 'UniformOutput' , false); 
    preceding_spike = cellfun(@(x) min(t-x) , preceding_spikes , 'UniformOutput' , false); 
    preceding_spike(cellfun('isempty', preceding_spike)) = {far_points(1)}; % maxing out to min time when I have no spike 
    %nonEmptyCells_preceding = preceding_spike(~cellfun('isempty', preceding_spike));
    preceding_spikes_array = cell2mat(preceding_spike); % This is tp(n)
    sigmap = std(preceding_spikes_array); 

    %Find following spikes 
    following_spikes = cellfun(@(x) x(find(x>t)) , thistrialtable , 'UniformOutput' , false); 
    following_spike = cellfun(@(x) max(t-x) , following_spikes , 'UniformOutput' , false); 
    following_spike(cellfun('isempty', following_spike)) = {far_points(2)};
    %nonEmptyCells_following = following_spike(~cellfun('isempty', following_spike));
    following_spikes_array = cell2mat(following_spike); % This is tp(n)
    sigmaf = std(following_spikes_array); 
    
    xp_m =1/ mean(t - preceding_spikes_array);
    xf_m = 1/mean(following_spikes_array -t);

    Num = sigmap*xf_m + sigmaf*xp_m; 
    Den = mean(following_spikes_array - preceding_spikes_array); 
    SPIKE(it) = Num/Den; 

end

% Plot spikes and synchrony measure
figure
hold on
%plot(tvect, SPIKE*1000, '-')
plot(time_points, SPIKE*10000, 'b-')
%
% plot(time_points, S_m_t*10000, 'r-')
hold on
for iu = 1 : length(thistrialtable)
    plot(thistrialtable{iu}, -iu*ones(size(thistrialtable{iu})), '|'); 
end
xlim([-0.1, 0.5])

%% 
multivariate_spiketrain = []
thistrialtable = {}

mintime = -1; 
maxtime = 1; 

dt = 0.001; %median(diff(sort(multivariate_spiketrain))); 

itrial = 6; 
for iu = 1 : size(all_mice_spikes_table, 1)
    multivariate_spiketrain = [multivariate_spiketrain; all_mice_spikes_table{iu, itrial}{1}]; 
    idx = intersect(find(all_mice_spikes_table{iu, itrial}{1}>mintime), find(all_mice_spikes_table{iu, itrial}{1}<maxtime)); 
    thistrialtable{iu} = all_mice_spikes_table{iu, itrial}{1}(idx); 
end
% Comparing pairwise and multivariate
mintime = -1; 
maxtime = 1; 
time_points = [mintime:dt:maxtime]; 

dt = 0.0001; 
% find cells with the most spikes
for iu = 1 : length(thistrialtable)
    numspikes(iu) = length(thistrialtable{iu}); 
end
[NS, I] = sort(numspikes); 
reduced_dataset = thistrialtable(I(end-20: end)); 
[D_Sm, S_m_t] = multivariate_average_spike_distance(reduced_dataset, max(timepoints) - min(time_points), time_points)

S = multivariateSPIKEDistance(reduced_dataset, time_points); 

figure
hold on
K = 10; 
plot(time_points, S_m_t*K , 'r-')
plot(time_points, S*K , 'b-')

hold on
for iu = 1 : length(reduced_dataset)
    plot(reduced_dataset{iu}, -iu*ones(size(reduced_dataset{iu})), '|'); 
end


%%
function K = sigmoidKernel(U, V)
    alpha = 0.01;  % scale
    c = 0;         % bias
    K = tanh(alpha * (U * V') + c);
end