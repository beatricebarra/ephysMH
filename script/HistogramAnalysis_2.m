% Reduced dataaset 
load('/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/24_11_18/reduced_dataset.mat')
analogFS = 1000; 
Sniffdt = 1/analogFS; 
spiketime_trains_data = spiketime_trains; 

spiketime_trains = {}; 

odornums =[5, 6, 7, 8]; 
odorlabels = {'ethyltiglate', 'ethylbutyrate', 'acetophenone', 'heptanone'}; 
%% Add sniff to dataset ( THIS DELETE ONCE YOU UPDATE DATASETS)
filename = '/Users/barrab01/Documents/PostDoc/project/Data/Monell_PCxEphys/Data/040118/24_11_18/040118__24_11_18.nwb'
nwbFile = nwbRead(filename);
trials_start_stop_time = [nwbFile.intervals.get('odor_trials').start_time.data(:), ...
                        nwbFile.intervals.get('odor_trials').stop_time.data(:)]; 

Odors = nwbFile.intervals.get('odor_trials').vectordata.get('pins').data(:); 
Concentration = nwbFile.intervals.get('odor_trials').vectordata.get('concentration').data(:); 

Sniff = nwbFile.acquisition.get('respiration_data').data(:);  
SniffTime = nwbFile.acquisition.get('respiration_data').timestamps(:);  
fSniff = lowPass(Sniff, round(1/Sniffdt), 10, 3); 
baselineSniff = lowPass(Sniff, round(1/Sniffdt),0.8, 2); 

plot_time = 2; 
for itrial = 1 : size(trials_start_stop_time, 1)
    itrial% find sniff snip 
    spiketime_trains{itrial} = spiketime_trains_data{itrial}; 
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
   
        % Plot
        plot(spiketime_trains{itrial}.SniffTime,filtered_sniff, 'k'  )
        plot(spiketime_trains{itrial}.SniffTime, logicSniff, 'm' )
        plot([min(spiketime_trains{itrial}.SniffTime), max(spiketime_trains{itrial}.SniffTime)], [th, th], 'r')
        plot(spiketime_trains{itrial}.inhalation_time, mean(spiketime_trains{itrial}.bSniff), 'ok', 'MarkerFaceColor','k')
         title(itrial)
        pause(0.011)
        clf
    catch
        spiketime_trains{itrial}.inhalation_idx = NaN; 
        spiketime_trains{itrial}.inhalation_time = NaN; 
        title('Trial without sniff')
        discarded_trials = [discarded_trials; itrial]; 
        plot(spiketime_trains{itrial}.SniffTime,filtered_sniff, 'k'  )
        plot(spiketime_trains{itrial}.SniffTime, logicSniff, 'm' )
        pause(0.011)        
        clf
    end
    
    
end

%% Extract trial numbers per condition
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


%% Find excitatory units PCx through PCA
% Plot PCA and kmeans analysis for PC and OB units

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
clustercols = {'r', 'g', 'b', 'y'}; 
[clusterPCCells, PCIDX] = cluster_units(PCMAT, clustercols); 
for i = unique(PCIDX)'
    PCClSize(i) = length(find(PCIDX ==i)); 
end
[U, exPCCL] = max(PCClSize); 

%% Working out histogram fitting
k = 3;
figure
for itrial = 1:length(spiketime_trains)
    
  
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


%% Generate gaussians for each trial and save parameters
% Also, plot averaged gausssians
k = 3; 
c = hot(12);
allGaussParams = []; 
conditionMat = []; 
histmat = []; 
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
            for iunit = goodPCunits(PCIDX == exPCCL) %% Piriform cortex  exPCCL
                idx = intersect(find(spiketime_trains{itrial}.units{iunit}> spiketime_trains{itrial}.inhalation_time+mint), ...
                    find(spiketime_trains{itrial}.units{iunit}< spiketime_trains{itrial}.inhalation_time + maxt)); 
                allunitsPC = [allunitsPC; spiketime_trains{itrial}.units{iunit}(idx)]; 
                
            end
        
            X = allunitsPC-(spiketime_trains{itrial}.inhalation_time); 
            
            if not(isempty(X)) % Sometimes there is no inhalation
                [mu, sigma] = GMM_clustering(X,k); 
                disp([mu, sigma])
                [mu, idx]=sort(mu); 
                sigma = sigma(idx); 
                disp([mu, sigma])

                allGaussParams = [allGaussParams; mu, sigma]; 
                conditionMat = [conditionMat; iodor, trials_numbers.(odorlabels{iodor}).concs(iconc)]; 
                histmat = [histmat; hist(X, 100)]; 
                MU = [MU; mu]; 
                SIGMA = [SIGMA; sigma]; 
                % Plotting 
                %figure(figPC)
                %subplot(1, length(trials_numbers.(odorlabels{iodor}).trial_idxs{iconc}), ii)
                %plot(spiketime_trains{itrial}.SniffTime - spiketime_trains{itrial}.inhalation_time, 100*spiketime_trains{itrial}.Sniff)
                %hold on
                %hist(X, 100)
                %hold on
                x = [0:0.001:1];
                % Plot over the existing figure, using black lines for the estimated pdfs.
                for i = 1 : k
            
                    y1 = gaussian1D(x, mu(i), sigma(i));
                    %plot(x, y1, '-','Color', c(iconc,:) , 'Linewidth', 0.5);
                end
               % ylim([0, 50])
            end
        end
        
        % Mean activation figure
        figure(figmean)
        hold on
        for i = 1 : 3
            y1 = gaussian1D(x, mean(MU(:,i)), mean(SIGMA(:, i)));
            plot(x, y1, '-','Color', c(iconc,:) , 'Linewidth', 2);
        end
        title([num2str(odorlabels{iodor}), '_' , num2str(trials_numbers.(odorlabels{iodor}).concs(iconc))])        
    end    
end

%% For each gaussian parameter, compute gaussian and extract features, plot them against concentration, different odorants
odor_symbols = ['o', 'x', 's', 'v']; 
x = [0:0.001:1]; % time vector
gfeats = []; 
for it = 1 : size(allGaussParams, 1)
    y1 = gaussian1D(x, allGaussParams(it, 1), allGaussParams(it, 4)); % this are parameters of the "first gaussian" = first spurt of spiking
    % Find features
    [peakV, peakloc] = max(y1); 
    peaktime = peakloc*mean(diff(x)); 
    spread = allGaussParams(it, 4); % for now I leave this as the pure std of the gaussian
    gfeats = [gfeats; peakV, peaktime, spread]; 

end

[COEFF, SCORE, LATENT, TSQUARED] = pca(gfeats); 
figure
subplot(1, 2, 1)
hold on
%plot(SCORE(:,1), SCORE(:,2), 'o')
for iodor = 1 : length(odorlabels)
    s = odor_symbols(iodor); 
    idx_odor = find(conditionMat(:,1)== iodor); 
    %plot(SCORE(idx,1), SCORE(idx,2), s)
    for iconc = 1 : length(trials_numbers.(odorlabels{iodor}).concs)
         idx_conc = find(conditionMat(:,2)== trials_numbers.(odorlabels{iodor}).concs(iconc));
        idx = intersect(idx_odor, idx_conc); 
        plot3(SCORE(idx,1), SCORE(idx,2), SCORE(idx,3), s, 'Color', c(iconc,:) , 'Linewidth', 2)
    end
end
subplot(1, 2, 2)
imagesc(COEFF)

%% Random forest regressor on concentration
idx_nonnegative = find(conditionMat(:,2)~=0); 
Y = log10(conditionMat(idx_nonnegative,2)); 
X = [gfeats(idx_nonnegative,1:2)]; 
%X = [SCORE(idx_nonnegative, 1:2)]; %this is worse
lab_odor = conditionMat(idx_nonnegative,1); 

figure;
hold on; 
for i = 1 : 5
    [rfModel, X_test, Y_test, Y_pred, mse] = emh_random_forest_regressor(X, Y, lab_odor, odor_symbols, '--'); 
    MSE_conc(i) = mse; 
end

%% Random forest regressor on intensity
I = zeros(size(conditionMat(:,1))) ; 
for iodor = 1:4
    idx_odor = (find(conditionMat(:,1)==iodor));
    concs = unique(conditionMat(idx_odor,2)); 
    for ic = 1:length(concs)
        idx_conc  = (find(conditionMat(:,2)==concs(ic)));
        idx = intersect(idx_odor, idx_conc); 
        if ic ==1
            I(idx)=  0; % no stim
        elseif ic ==2 || ic ==3
            I(idx)=  1;  %very low
        elseif ic ==4 
           I(idx)=  2;  % low
        elseif ic ==5 
           I(idx)=  3;  % medium
        elseif ic ==6 
           I(idx)=  4; % high
        else 
           I(idx)=  5;  % very high
        end
    end
end

idx_nonnegative = find(conditionMat(:,2)~=0); 
Y = I (idx_nonnegative); 
X = [gfeats(idx_nonnegative,1:2)]; 
%X = [SCORE(idx_nonnegative, 1:2)]; %this is worse
lab_odor = conditionMat(idx_nonnegative,1); 
figure;
hold on; 
for i = 1 : 5
    [rfModel, X_test, Y_test, Y_pred, mse] = emh_random_forest_regressor(X, Y, lab_odor, odor_symbols); 
    MSE_conc(i) = mse; 
end

%% Reducing to first gaussian is not enough to perform well on intensity categories
% I am gonna try to feed the histogram bins and see what happens

%% Histograms predictions of concentration
idx_nonnegative = find(conditionMat(:,2)~=0); 
Y = log10(conditionMat(idx_nonnegative,2)); 
X = [histmat(idx_nonnegative,:)]; 
%X = [SCORE(idx_nonnegative, 1:2)]; %this is worse
lab_odor = conditionMat(idx_nonnegative,1); 
figure;
hold on; 
for i = 1 : 5
    [rfModel, X_test, Y_test, Y_pred, mse] = emh_random_forest_regressor(X, Y, lab_odor, odor_symbols); 
    MSE_conc(i) = mse; 
end
 
 %% Histograms predictions of intensity "classes"
idx_nonnegative = find(conditionMat(:,2)~=0); 
Y = I(idx_nonnegative); 
X = [histmat(idx_nonnegative,:)]; 
%X = [SCORE(idx_nonnegative, 1:2)]; %this is worse
lab_odor = conditionMat(idx_nonnegative,1); 
figure;
hold on; 
for i = 1 : 5
    [rfModel, X_test, Y_test, Y_pred, mse] = emh_random_forest_regressor(X, Y, lab_odor, odor_symbols); 
    MSE_Int(i) = mse; 
end
%% Histograms/gaussian predictions of intensity as per perceived intensity curves

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
for ii = 1 : length(conditionMat(:,2))
    % convert concentration ot intensity
    io = conditionMat(ii, 1); 
    disp(log10(conditionMat(ii, 2)*parammat(io, 4)))
    realI(ii) = psychocurve(log10(conditionMat(ii, 2)*parammat(io, 4)), parammat(io, 1), parammat(io, 2), parammat(io, 3)); 
    plot(log10(conditionMat(ii, 2)*parammat(io, 4)), realI(ii), 'o', 'MarkerFaceColor', int_curves_cols(io)) 

end
Y = realI(idx_nonnegative)'; 
X = [histmat(idx_nonnegative,:)]; 
%X = [gfeats(idx_nonnegative,:)];
X = allGaussParams(idx_nonnegative,[4:6]); 
lab_odor = conditionMat(idx_nonnegative,1); 
figure;
hold on; 
MSE_realInt = []; 
for i = 1 : 5
    [rfModel, X_test, Y_test, Y_pred, mse] = emh_random_forest_regressor(X, Y, lab_odor, odor_symbols, '-'); 
    MSE_realInt(i) = mse; 
end




%% Not great, when you compute mape you realize that the error is huge with respect to the scale

% Visualize averages
figure
hold on
x = linspace(0, 0.8, 100); 
for iodor = 1 : 4
   
    idx_odor = (find(conditionMat(:,1)==iodor));
    for ii = 2:4
        idx_int = find(I == ii); 
        idx = intersect(idx_odor, idx_int); 
        P = mean(allGaussParams(idx,:)); 
        subplot(1, 3, ii-1)
        hold on
        y1 = gaussian1D(x, P(1), P(4));
        plot(x, y1, int_curves_cols(iodor), 'Linewidth', 2)
        ylim([0, 24])
        xlim([0, 0.2])
  
    end
end


%%
figure
hold on
for io = 1 : 4
x = logspace(-1, 6, 100); 
x = log10(x); 
y = psychocurve(x, parammat(io, 1),  parammat(io, 2), parammat(io, 3)); 

plot(x, y, int_curves_cols(io),  'Linewidth', 2) 
%set(gca,'Xscale', 'log')
end
function [y] = psychocurve(c, a, k, n)
    disp([a,k,n])
    for i = 1 : length(c)
        y(i) = a/(1 + exp(-k*(c(i)-n))); 
    end
    %y = a./(1 + exp(-k.*(c-n))); 
end