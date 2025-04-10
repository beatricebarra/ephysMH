function [OB_units_idx, PC_units_idx] = units_quality_check(nwbFile)
    
location_drop = 'not in';  % Example: location filter for specific areas
    ks_label = 'good';  % Example: Kilosort quality
    SI_label = 'good';  % Example: SI quality
    cutoff_hz = 0.2;  % Example: firing rate cutoff
    cutoff_amp = 10 ; % Example: amplitude cutoff
    max_above = 400 ; % Example: max depth for 'PC:' locations
    % Bea modified to 
    cutoff_amp = 3 ; % Example: amplitude cutoff
    max_above = 400 ; % Example: max depth for 'PC:' locations
    
    % Compute depth
    temp = cellfun(@(x) split(x, ':'), nwbFile.units.vectordata.get('location').data(:), 'UniformOutput', false); 
    depth = cellfun(@(x) str2num(x{2}), temp); 
    
    % General quality check (same for PC and OB)
    quality_units = (cellfun(@(x) ~contains( x , location_drop ), nwbFile.units.vectordata.get('location').data(:)) & ...
        cellfun(@(x) contains( x , ks_label ), nwbFile.units.vectordata.get('quality').data(:)) & ...
        cellfun(@(x) contains( x , SI_label ), nwbFile.units.vectordata.get('SI_quality').data(:)) &...
        nwbFile.units.vectordata.get('firing_rate').data(:)>cutoff_hz & ...
        (nwbFile.units.vectordata.get('amplitude').data(:)>cutoff_amp)'); 
    
    PC_units = (cellfun(@(x) contains( x , 'PC' ), nwbFile.units.vectordata.get('location').data(:)) & ...
        depth<max_above & ...
        quality_units);
    PC_units_idx = find(PC_units); 
    OB_units = (cellfun(@(x) contains( x , 'OB' ), nwbFile.units.vectordata.get('location').data(:)) & ...
        quality_units); 
    OB_units_idx = find(OB_units); 
    %disp([sum(PC_units), sum(OB_units)])
end