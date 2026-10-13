function [amp, baseline_condition] = ComputeAmplitude( ts, para, varargin)
stats = [];
epoch_time = [];
ev_epochs = [];
comp_data = [];
sd = [];
amp = [];
nbpermuts = 10000;


if ~isempty( varargin)
    amp = abs( varargin{1});
    nbsamples = size( amp, 2);
    comp_data = varargin{1};
end


events = [];
if isfield( para, 'events')
    events = para.events;
end

epoch = [];
if isfield( para, 'epoch')
    epoch = para.epoch;
end

smooth_win = [];
if isfield( para, 'smooth_win')
    smooth_win = para.smooth_win;
end

if isfield( para, 'nbpermuts')    
    nbpermuts = para.nbpermuts;
end

freq = para.freq;
nbfreqs = length( freq);
bandwidth = para.bandwidth;
srate = para.srate;

if ~isfield( para, 'measure')
    measure = 'amplitude';
    
else
    measure = para.measure;
    switch measure
        case {'amplitude', 'power'}
            % do nothing
        otherwise
            error( 'Unknown measurement.');
    end
end

if ~isfield( para, 'trl')
    para.trl = 'mean';
end

if isempty( amp)
    if ~isvector( ts)
        error( 'Time series must be a vector.');
    end

    nbsamples = length( ts);
    ts = reshape( ts, [1, nbsamples]);    
    comp_data = zeros( nbfreqs, nbsamples);
    for freq_idx = 1 : nbfreqs  
        bpfreq = freq( freq_idx) + [-1, 1]*bandwidth( freq_idx)/2;
        fprintf( 'Bandpass filtering %d/%d: %.1f - %.1f Hz               \r', freq_idx, nbfreqs, bpfreq(1), bpfreq(2));
        if bpfreq(1) < eps
            filt_ts = ft_preproc_lowpassfilter( ts, srate, bpfreq(2), [], 'fir');
        else
            filt_ts = ft_preproc_bandpassfilter( ts, srate, bpfreq, [], 'fir');
        end
        
        comp_data( freq_idx, :) = hilbert( filt_ts);
    end
    
    % amplitude
    amp = abs( comp_data);    
    fprintf( 'Bandpass filtering                                       Done.\n');
end
% amplitude -> power
if strcmpi( measure, 'power')
    amp = amp .^2;
end   

% temporal smoothing
if ~isempty( smooth_win)
    smooth_kernel = round( smooth_win * srate);
    if mod( smooth_kernel, 2) == 0
        smooth_kernel = smooth_kernel + 1;
    end
    amp = G_Smooth( amp, 2, smooth_kernel);
end
    
if isempty( events)
    return;
end

if isempty( epoch) || length( epoch) ~= 2
    error( 'Epoch must be given as [lo, hi] in seconds.');
end

if ~isfield( para, 'baseline_condition')
    baseline_condition = [];
else
    baseline_condition = para.baseline_condition;
    if ~isempty( baseline_condition)
        if ~ismember( baseline_condition, 1:length(events))
            error( 'Baseline condition is out of range.');
        end
    end
end

end