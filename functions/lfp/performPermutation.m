function [stats] =performPermutation(amp, resp, nbevs, para, ts, varargin)
    events = [];
    if isfield( para, 'events')
        events = para.events;
    end
    
    epoch = [];
    if isfield( para, 'epoch')
        epoch = para.epoch;
    end
    nbsamples = length( ts);
    srate = para.srate;

    if ~isempty( varargin)
        amp = abs( varargin{1});
        nbsamples = size( amp, 2);
        comp_data = varargin{1};
    end
    freq = para.freq;
    nbfreqs = length( freq);
    if isfield( para, 'nbpermuts')    
        nbpermuts = para.nbpermuts;
    end
    
    sd = cell( nbevs, 1);
    
    h = waitbar(0, 'Please wait...');
    for ev_idx = 1 : nbevs 
        % Update waitbar
        waitbar(ev_idx/nbevs, h, sprintf('Progress: %d%%', ev_idx));
        if ~isempty( events{ ev_idx})    
    
            current_event = cast( events{ ev_idx}, 'like', nbsamples);
            permut_avg = zeros( nbfreqs, nbpermuts);
            for pidx = 1 : nbpermuts
                %fprintf( 'Permutation %d/%d               \r', pidx, nbpermuts);
    %             shift_amount = randi( [int32( 0.1*srate), nbsamples], 1);
                shift_amount = randi( [cast( 0.1*srate, 'like', nbsamples), nbsamples], 1);
                permut_onsets = mod( current_event + shift_amount, nbsamples);
                % make sure there' enough sample at the beginning and end
                while any( permut_onsets(:) < -1*epoch(1)*srate + 1) || any( permut_onsets(:) > nbsamples - epoch(2)*srate - 1)
                    %shift_amount = randi( [int32( 0.1*srate), nbsamples], 1);
                    shift_amount = randi( [cast( 0.1*srate, 'like', nbsamples), nbsamples], 1);
                    permut_onsets = mod( current_event + shift_amount, nbsamples);
                end
                switch para.trl
                    case 'mean'
                        permut_avg( :, pidx) = mean( amp( :, permut_onsets), 2);
                    otherwise
                        permut_avg( :, pidx) = median( amp( :, permut_onsets), 2);
                end
            end
    
            sd{ ev_idx} = std( permut_avg, [], 2);
            Z{ ev_idx} = bsxfun( @rdivide, resp{ ev_idx}, sd{ ev_idx});
        
        else
            Z{ ev_idx} = [];
        end
    end % event 
    fprintf( 'Permutation                           Done.\n');
    
    % stats.resp_percentchange = resp_percentchange;
    stats.Z = Z;
    % stats.resp = resp;
end