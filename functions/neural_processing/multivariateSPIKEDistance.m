function S = multivariateSPIKEDistance(spikes, t)
% Computes the true multivariate SPIKE-distance
% spikes : cell array of spike trains, e.g. {spk1, spk2, ...}
% t      : vector of time points at which to compute the distance

N = numel(spikes);
t = t(:)'; % ensure row vector
T_start = min(t); 
T_end = max(t); 
T = length(t);
S = zeros(1, T); % preallocate

% Handle trivial case: no spike trains
if N == 0
    S = zeros(1, T);
    return;
end

% Precompute previous and next spike times for all neurons
t_prev = NaN(N, T);
t_next = NaN(N, T);
isi = NaN(N, T); % ISI at each time
countnorm = 0; 
for n = 1:N
    spk = sort(spikes{n}(:)');
    
    if isempty(spk)
        % No spikes -> set previous to -Inf, next to +Inf
        t_prev(n, :) = -Inf;
        t_next(n, :) = Inf;
        isi(n, :) = Inf;
    elseif numel(spk) == 1
        % Only one spike -> set prev = next = that spike
        t_prev(n, :) = spk;
        t_next(n, :) = spk;
        isi(n, :) = Inf; % ISI undefined -> infinite
    else
        % Normal case
        for it = 1 : length(t)
            ti = t(it); 
            %t_prev(n,it) = interp1(spk, spk, t, 'previous', 'extrap');
            %t_next(n,:) = interp1(spk, spk, t, 'next', 'extrap');
            % Previous spike or T_start
            if isempty(spk(spk <= ti))
                t_prev(n, it) = T_start; 
            else
                t_prev(n, it) = max(spk(spk <= ti));
            end
            if isempty(spk(spk > ti))
                t_next(n, it) = T_end; 
            else
                t_next(n, it) = min(spk(spk > ti)); 
            end
        end
        isi(n,:) = t_next(n,:) - t_prev(n,:);
    end
end

% Compute distance profile
for i = 1:T
    dt_prev = abs(t(i) - t_prev(:,i));
    dt_next = abs(t_next(:,i) - t(i));
    
    % Local spike time difference for each neuron
    s_local = min(dt_prev, dt_next) ./ isi(:,i);
    
    % Handle division by Inf (empty or singleton spike trains)
    s_local(isinf(s_local)) = 0;
    
    % The true multivariate SPIKE distance is the average across neurons
    S(i) = nansum(s_local) / N;
end
end
