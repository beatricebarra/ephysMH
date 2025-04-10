function [FRmat, peak_latency ]= extract_peak_latency(spikes_table, sig, mint, maxt)
tvect = linspace(mint, maxt, (maxt-mint)/sig); 
tvectbaseline = linspace(mint-10*sig, mint-sig, (mint-sig-(mint-10*sig))/sig); 
%tvectall = linspace(mint-10*sig, maxt, (maxt-(mint-10*sig))/sig); 
chronuxfolder = '/Users/barrab01/Documents/Repos/chronux_2_12'
for icondition = 1 : size(spikes_table, 2)
       
        FRmat_cond = []; 
        sorted_responses = []; 
        for iunit = 1 : size(spikes_table, 1)
            [FR,t, err] = compute_FR_from_raster(spikes_table{iunit, icondition} , sig, mint, maxt, tvect ); 
            %[baselineFR,t, err] = compute_FR_from_raster(spikes_table{iunit, icondition} , sig, mint-10*sig, mint-sig, tvectbaseline ); 
            FRmat_cond= [FRmat_cond; FR]; % append to create matrix with all responses
            %[allFR,t, err] = compute_FR_from_raster(spikes_table{iunit, icondition} , sig, mint-10*sig, maxt, tvectall ); 
            % with chronux peak finder
            %[xmax] = findpeaks(FR(2:end), 1)
            %xmax.loc = xmax.loc + 1; 
            %[maxidx]= find_true_peak(FR, xmax); 

            % with normal peak finder
            %rmpath(genpath(chronuxfolder));
            [peaks, locs] = findpeaks_Matlab(FR, 'MinPeakHeight', 5, 'MinPeakProminence', 3); %mean(baselineFR)+std(baselineFR)
            
            xmax.loc = locs; 
            [maxidx]= find_true_peak(FR, xmax); 
%             plot(FR)
%             hold on
%             plot(locs, peaks, 'or')
%             if isnan(maxidx)
%                title('NaN')
%             else    
%                 plot(maxidx, FR(maxidx), 'ob')
%             end
%             pause
%             clf
           
            %addpath(genpath(chronuxfolder));
            

            %[maxv, maxloc] = max(FR(xmax.loc)); 
            % Test plot
            %hold on
            %plot(tvectall, allFR)
            %plot(tvect, FR)
            %plot([min(tvect), max(tvect)], [mean(baselineFR)+std(baselineFR), mean(baselineFR)+std(baselineFR)], 'r--')
            %maxidx
            %if isnan(maxidx)
            %    title('NaN')
            %else    
            %    plot(tvect(maxidx), FR(maxidx), 'o')
            %end
            
            peak_latency(iunit, icondition) = (maxidx-1)*sig; 
            
            
        end
        disp('Checking stuff in extract peak latency ')
        
        FRmat{icondition} = FRmat_cond; 
    end
end