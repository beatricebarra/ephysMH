function [D_Sm, S_m_t] = multivariate_average_spike_distance(spike_trains, T, time_points)
% Compute the multivariate SPIKE-distance for a set of spike trains.
%
% Input:
%   spike_trains - Cell array of spike time vectors, one per neuron
%   T            - Total duration of recording
%
% Output:
%   D_Sm         - Multivariate SPIKE-distance

num_neurons = length(spike_trains);

if num_neurons == 0
    D_Sm = 0;
    return;
end

% Instantaneous SPIKE-distance profile function
function Smt = instantaneous_spike_distance(t)
    distances = [];
    for i = 1:num_neurons
        for j = i+1:num_neurons
            % Find nearest spike to t for neuron i
            [~, idx_i] = min(abs(spike_trains{i} - t));
            nearest_i = spike_trains{i}(idx_i);
    
            % Find nearest spike to t for neuron j
            [~, idx_j] = min(abs(spike_trains{j} - t));
            nearest_j = spike_trains{j}(idx_j);
    
            % Absolute difference
            distances = [distances, abs(nearest_i - nearest_j)];
        end
    end
    if isempty(distances)
        Smt = 0;
    else
        Smt = mean(distances);
    end
end

% Integration over time
 % number of evaluation points
S_m_t = arrayfun(@instantaneous_spike_distance, time_points);

% Trapezoidal integration
D_Sm = trapz(time_points, S_m_t) / T;

end