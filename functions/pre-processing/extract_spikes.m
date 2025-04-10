function [spikes_table, count_table, condition_table] = extract_spikes(nwbFile, spikes, Odors, Concentration, Prex, Postx, FVO,  PC_units_idx, all_odors)
    Nodors = length(unique(Odors)); 
    Nconcall = []; 
    for iodor = unique(Odors)'
        idx = find(Odors == iodor); 
        Nconcall= [Nconcall, length(unique(Concentration(idx)))]; 
    end

    Nconcs = mean(Nconcall); 
    disp('Printing concentration numbers')
    [Nconcall, Nconcs]
    mywindow = 0.5; 
    spikes_table = cell(length(PC_units_idx), Nodors*Nconcs);
    count_table = zeros(length(PC_units_idx), Nodors*Nconcs); 
    condition_table = zeros( Nodors*Nconcs, 4); 
    for itrial  = 1 : length(Odors)
        
        % Find odor and concentration indexes to pull together trials with the
        % same conditions
        iodor = find(all_odors == Odors(itrial)); 
        idx_odor = find(Odors ==iodor+4); 
        all_conc_thisodor = unique(Concentration(idx_odor));
        iconc = find(all_conc_thisodor == Concentration(itrial)); 
        conc_value = Concentration(itrial); 
        condition_idx = (iodor-1)*length(all_conc_thisodor) + iconc; 
        condition_table( condition_idx, 4) = condition_table( condition_idx, 4) +1; 
        condition_table( condition_idx, 1:3) = [iodor, iconc, conc_value]; 
    
        all_prex_before_fvo_idx = find(Prex<FVO(itrial));
        all_prex_after_fvo_idx = find(Prex>FVO(itrial));
        all_postx_before_fvo_idx=  find(Postx<FVO(itrial));
        
        % Start and end of baseline
        start_bas = Prex(all_prex_before_fvo_idx(end-2));  
        end_bas = start_bas + mywindow; Postx(all_postx_before_fvo_idx(end)); 
        
        % Start and end of inh
        start_inh = Prex(all_prex_after_fvo_idx(1)); 
        end_inh =start_inh + mywindow; 
        
        % Loop through units and find spikes in baseline and first sniff
        iiunit = 0; 
        for iunit = PC_units_idx'
            iiunit = iiunit+1; 
            last_spike = nwbFile.units.spike_times_index.data(iunit); 
            if iunit >1
                first_spike = nwbFile.units.spike_times_index.data(iunit-1)+1; 
            else
                first_spike = 1; 
            end
    
            allspikes = spikes(first_spike:last_spike); 
            allspikesintrial = allspikes(intersect(find(allspikes>start_bas), find(allspikes<start_inh+2))); 
            if isempty(spikes_table{iiunit, condition_idx})
                spikes_table{iiunit, condition_idx} = {allspikesintrial-start_inh}; 
                count_table(iiunit, condition_idx) = count_table(iiunit, condition_idx) +1; 
            else
                spikes_table{iiunit, condition_idx}{count_table(iiunit, condition_idx) +1} = allspikesintrial-start_inh; 
                count_table(iiunit, condition_idx) = count_table(iiunit, condition_idx) +1;
            end
        end
        percent_complete = (itrial / length(Odors)) * 100;
        disp([' ----   ', num2str(percent_complete), ' % of Extracting spikes   ----'])
    end
end