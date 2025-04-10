function [Synch] = computeSynchMeasure(spikes, time_points, maxW)
    
    for it = 1 : length(time_points) % for each time instant
        closest_spike_array = [];     
        t = time_points(it ); % Compute time 
        
        [value_max, idx_closest_spike] = cellfun(@(x) min(abs(x-t)) , spikes , 'UniformOutput' , false); 
        
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
        closest_spike_array = cell2mat(new_array); 
        %hist(closest_spike_array)
        %title([num2str(t), '_', num2str(1- mean(closest_spike_array)/(maxW))])
        
        Synch(it) = 1- sum(closest_spike_array./(maxW))/length(spikes); 
    end