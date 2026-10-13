function [resp, resp_percentchange, nbevs] = correctBaseline(amp, baseline_condition, para, varargin)

events = [];
if isfield( para, 'events')
    events = para.events;
end

epoch = [];
if isfield( para, 'epoch')
    epoch = para.epoch;
end

freq = para.freq;
nbfreqs = length( freq);
bandwidth = para.bandwidth;
srate = para.srate;


% -------
seg = floor( epoch * srate);
nbevs = length( events);
resp = cell( nbevs, 1);
resp_percentchange = cell( nbevs, 1);
Z = cell( nbevs, 1);
ev_epochs = cell( nbevs, 1);
for ev_idx = 1 : nbevs
    current_event = events{ ev_idx};
    if isempty( current_event)
        ev_epochs{ ev_idx} = [];
        resp{ ev_idx} = [];
        
    else
        % freq, samples x epochs
        ev_epochs{ ev_idx} = G_EpochTS( amp, current_event, seg);
        % average across trials
        switch para.trl
            case 'mean'
                resp{ ev_idx} = mean( ev_epochs{ev_idx}, 3);
            case 'median'
                resp{ ev_idx} = median( ev_epochs{ev_idx}, 3);
            otherwise
                error( 'error');
        end  
    end
end % event

% baseline correction
epoch_time = epoch(1) : 1/srate : epoch(2);
baselineloc = [];
if isfield( para, 'baseline')
    baseline = para.baseline;
    baselineloc = G_RangeLoc( epoch_time, baseline);
end

if ~isempty( baseline_condition)
    fprintf( 'Using condition %d as a common baseline\n\n', baseline_condition);
    if ~isempty( events{ baseline_condition})
        avg = mean( resp{ baseline_condition}( :, baselineloc, :), 2);
        for ev_idx = 1 : nbevs
            resp{ ev_idx} = bsxfun( @minus, resp{ ev_idx}, avg);
            resp_percentchange{ ev_idx} = 100* bsxfun( @rdivide, resp{ ev_idx}, avg);
        end % event
        
    else
        error( 'Baseline condition has no events.');
    end

else
    fprintf( 'Perform baseline correction for each condition separately.\n\n');
    for ev_idx = 1 : nbevs
        if ~isempty( events{ ev_idx})
            avg = mean( resp{ ev_idx}( :, baselineloc, :), 2);
            resp{ ev_idx} = bsxfun( @minus, resp{ ev_idx}, avg);
            resp_percentchange{ ev_idx} = 100* bsxfun( @rdivide, resp{ ev_idx}, avg);
        else
            resp_percentchange{ ev_idx} = [];
        end
        
    end % event
end

end