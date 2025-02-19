function [sortedFRMAT,Isort, sortedv]  = sort_cells_by_peak(FRMAT, tvect, Nstd)

%THIS HAS A PROBLEM
    crossth = []; 
    trigavgmat = FRMAT; 
    mint = tvect(1); 
    dt = mean(diff(tvect)); 
    for  icell = 1 : size(trigavgmat, 1)
        % is there a peak for that cell 
        [peaks, peaklocs] = findpeaks(trigavgmat(icell,:)); 
        idx_real_peaks = intersect(find(peaklocs> (0-mint)/dt), find(peaklocs<0.5/dt)); 
        idx_real_peaks = intersect(idx_real_peaks, find(peaks>mean((trigavgmat(icell, :))) + Nstd*std((trigavgmat(icell, :))))); 
        if ~isempty(idx_real_peaks)
            realpeaks = peaklocs(idx_real_peaks); 
            crossth(icell) = realpeaks(1); 
        else
            crossth(icell) = NaN; 
        end
    end    
    idxzeros = find(crossth == 0); 
    crossth(idxzeros) = NaN; 
    

    [sortedv, Isort] = sort(crossth);
    sortedFRMAT = []; 
    iicell = 0; 
    for icell = Isort
        iicell = iicell +1; 
        sortedFRMAT(iicell,:) = FRMAT(icell,:); 
    
    end