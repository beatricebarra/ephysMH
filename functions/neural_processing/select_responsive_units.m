function [selunits] = select_responsive_units(all_mice_spikes_table, all_mice_condition_table, iunits)
    
    condition_table = all_mice_condition_table(1:32,:); 
    odors = unique(condition_table(:,1)); 
    selunits = {}; 
    for iodor = 1 : length(unique(condition_table(:,1)))
        selunits{iodor} = []; 
    end
    for iu = iunits'
        %disp(['Unit ', num2str(iu)])
        for iodor = 1 : length(unique(condition_table(:,1)))
            idxO = find(condition_table(:,1) == odors(iodor)); 
            concs = unique(condition_table(idxO,3));
            concs = concs(concs>0); 
            for iconc = 1: length(concs)
                idxC = find(condition_table(:,3) == concs(iconc)); 
                idxT = intersect(idxO, idxC); 
                
                % Compute whether unit is responsive in this condition
                [response_index, pvalue] = find_responsive_units( all_mice_spikes_table{iu, idxT}, [0, 0.15], [-0.2, -0.05]); 
                responsiveness(iu, iodor, iconc) = pvalue;
    
            end
            % Find if there is a response for any of the concentrations on
            % this odorant
            if any(responsiveness(iu, iodor, :)<0.05) % save unit as viable only if responsive
                selunits{iodor} = [selunits{iodor}, iu] ; 
            end
        end
    end
end