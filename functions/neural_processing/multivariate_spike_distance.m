function S_t = multivariate_spike_distance(t, spike_trains)

% spike_trains: cell array where spike_trains{i} contains spike times of neuron i
% t: vector of time points at which to evaluate the distance

N = numel(spike_trains); % number of neurons
M = numel(t);            % number of time points
D = zeros(1, M);         % output distance

% Precompute previous and next spike times for each neuron
t_prev = NaN(N, M);
t_next = NaN(N, M);

for n = 1:N
    spk = spike_trains{n}
%     t_prev(n,:) = interp1(spk, spk, t, 'previous', 'extrap');
%     t_next(n,:) = interp1(spk, spk, t, 'next', 'extrap');
end

% Compute instantaneous spike distance at each time point
for m = 1:M
    % Get t_prev and t_next for all neurons at time t(m)
    tp = t_prev(:, m);
    tn = t_next(:, m);
    
    % Compute ISI for each neuron
    isi = tn - tp;
    
    % Avoid division by zero
    isi(isi == 0) = eps;
    
    % Compute distances from t(m) to previous and next spike for all neurons
    dtp = abs(t(m) - tp);
    dtn = abs(tn - t(m));
    
    % Local spike timing difference for each neuron
    S_local = min(dtp, dtn) ./ isi;
    
    % Average over all neurons
    D(m) = mean(S_local);
end

end