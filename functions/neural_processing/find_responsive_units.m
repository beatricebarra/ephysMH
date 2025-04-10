function [response_index, p] = find_responsive_units(spike_trains, inh, baseline)
% spike_train = cell array with a series of arrays contining spike times all spike times
% inh = arra with times for start and end inhalation
% baseline = array with times for start and end baseline
    start_inh = inh(1); 
    end_inh = inh(2); 
    inhwindow = end_inh - start_inh ; 
    start_bas = baseline(1); 
    end_bas = baseline(2); 
    baswindow = end_bas- start_bas; 
    stimulus_rates = []; 
    baseline_rates = []; 
    for itrain = 1 : length(spike_trains)
        % First sniff 
        stimulus_rates = [stimulus_rates, length(find( spike_trains{itrain}> start_inh & spike_trains{itrain}< end_inh))/ inhwindow];  
        % Baseline
        baseline_rates =[ baseline_rates, length(find( spike_trains{itrain}> start_bas & spike_trains{itrain}< end_bas))/ baswindow]; 
    end
    
    % Compute response index
    labels = [repmat([0], 1, length(baseline_rates)) ,  repmat([1], 1,length(stimulus_rates))];  
    scores = [baseline_rates , stimulus_rates];    
    %[X,Y,T,AUC] = perfcurve(labels,scores, 1); 
    [X, Y, AUC, p] = ranksumROC(labels,scores); 
    %plot(X, Y)
    %hold on
    response_index = 2 * (AUC - 0.5);
end