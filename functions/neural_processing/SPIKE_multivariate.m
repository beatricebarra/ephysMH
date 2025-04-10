function D = SPIKE_multivariate(spikes, T_start, T_end, dt)
% True multivariate SPIKE-distance
% Robust to empty spike trains or sparse spiking

N = length(spikes);
t_sample = T_start:dt:T_end;
D = zeros(size(t_sample));

dt = mean(diff(t_sample)); 
total = floor((max(t_sample)-min(t_sample))/dt); 
checkpoints = floor(linspace(0.01, 1, 100) * total);  % 10%, 20%, ..., 100%

for i = 1:length(t_sample)
    if any(i== checkpoints)
        fprintf('Progress: %.0f%% complete\n', (i / total) * 100);
        pause(0.01)
    end
    ti = t_sample(i);
    prev_spike = zeros(N,1);
    next_spike = zeros(N,1);
    isi = zeros(N,1);
    S_local = zeros(N,1);
    
    for n = 1:N
        s = spikes{n};
        
        if isempty(s)
            % Fully empty train -> maximal dissimilarity
            prev_spike(n) = T_start;
            next_spike(n) = T_end;
            isi(n) = T_end - T_start;
            S_local(n) = 1;  % Maximal dissimilarity
            continue
        end
        
        % Previous spike or T_start
        prev_spike(n) = max(s(s <= ti), T_start);
        
        % Next spike or T_end
        next_spike(n) = min(s(s > ti), T_end);
        
        % Local ISI
        isi(n) = next_spike(n) - prev_spike(n);
        
        % Avoid division by zero
        if isi(n) == 0
            S_local(n) = 0;
        else
            x_prev = ti - prev_spike(n);
            x_next = next_spike(n) - ti;
            S_local(n) = abs(x_next - x_prev) / isi(n);
        end
    end
    
    % Multivariate profile = mean of local profiles
    D(i) = mean(S_local);
end

end