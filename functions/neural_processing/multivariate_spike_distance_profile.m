function S = multivariate_spike_distance_profile(tVec, spikes)
% Compute the multivariate SPIKE-distance over a time vector
% tVec   - vector of times to evaluate the distance
% spikes - cell array of spike trains

S = zeros(size(tVec));
dt = mean(diff(tVec)); 
total = floor((max(tVec)-min(tVec))/dt); 
checkpoints = floor(linspace(0.01, 1, 100) * total);  % 10%, 20%, ..., 100%


for k = 1:length(tVec)
    if any(k == checkpoints)
        fprintf('Progress: %.0f%% complete\n', (k / total) * 100);
        pause(0.01)
    end
    S(k) = multivariate_spike_distance(tVec(k), spikes);
end

end