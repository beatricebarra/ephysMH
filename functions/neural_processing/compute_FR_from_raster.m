function [FR,t, err] = compute_FR_from_raster(Rast, sig, mint, maxt, tvect )
             
    % CReate structure for FR computation
    mystruct = struct(); 
    ii = 0; 
    for i = 1 : length(Rast)
        if isempty(Rast{i})
            disp('empty')
        else
            ii = ii +1 ; 
            mystruct(ii).FR = Rast{i}; 
        end
        
    end

    % Compute psth for this specific cell
    if isfield(mystruct, 'FR')
    
        [FR,t, err] = psth(mystruct, sig,'n', [mint, maxt], 2, tvect); % 2 means bootstrap
    else 
        FR = zeros(1, length(tvect)); 
        err = zeros(1, length(tvect)); 
        t = tvect; 
    end
end