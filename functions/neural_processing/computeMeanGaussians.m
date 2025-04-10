function [MeanGaussians , x, UpperboundGaussians, LowerboundGaussians]= computeMeanGaussians(spiketime_trains, trials_numbers, goodunits, mint, maxt,  k)

allGaussParams = []; 
conditionMat = []; 
histmat = []; 
odorlabels = fieldnames(trials_numbers);
x = [mint:0.001:maxt];
for iodor = 1: length(odorlabels)
    for iconc = 1 : length(trials_numbers.(odorlabels{iodor}).concs)
        
        % Initialize structure
        MeanGaussians{iodor, iconc} = cell(k,1);
        UpperboundGaussians{iodor, iconc} = cell(k,1);
        LowerboundGaussians{iodor, iconc} = cell(k,1);
        for igauss = 1 : k
            MeanGaussians{iodor, iconc}{igauss} = []; 
            UpperboundGaussians{iodor, iconc}{igauss} = []; 
            LowerboundGaussians{iodor, iconc}{igauss} = []; 
        end

        MU = []; 
        SIGMA = []; 
        for ii = 1 : length(trials_numbers.(odorlabels{iodor}).trial_idxs{iconc})
            allunitsPC = [];
            itrial = trials_numbers.(odorlabels{iodor}).trial_idxs{iconc}(ii); 
            goodunits
            for iunit = goodunits %% Piriform cortex  exPCCL
                idx = intersect(find(spiketime_trains{itrial}.units{iunit}> spiketime_trains{itrial}.inhalation_time+mint), ...
                    find(spiketime_trains{itrial}.units{iunit}< spiketime_trains{itrial}.inhalation_time + maxt)); 
                allunitsPC = [allunitsPC; spiketime_trains{itrial}.units{iunit}(idx)]; 
                
            end
        
            X = allunitsPC-(spiketime_trains{itrial}.inhalation_time); 
            
            if not(isempty(X)) % Sometimes there is no inhalation
                [mu, sigma] = GMM_clustering(X,k); 
                [mu, idx]=sort(mu); 
                sigma = sigma(idx);  
                MU = [MU; mu]; 
                SIGMA = [SIGMA; sigma]; 
%                 for i = 1 : k
%             
%                     y1 = gaussian1D(x, mu(i), sigma(i));
%                     if i ==1
%                         plot(x, y1, '-', 'Linewidth', 0.5);
%                         hold on
%                         
%                     end
%                     MeanGaussians{iodor, iconc}{i}(ii, :)= y1; 
%                 end
%                
%                ylim([0, 50])
            end
        end
        for i = 1 : k
            y1 = gaussian1D(x, mean(MU(:, i)), mean(SIGMA(:, i)));
            MeanGaussians{iodor, iconc}{i}= y1; 

            noise_mu = std(MU(:, i));%sqrt((1/(size(MU, 1)-1))*sum( abs(MU(:,i) - mean(MU(:, i))) )); 
            noise_sigma = std(SIGMA(:, i));%sqrt((1/(size(SIGMA, 1)-1))*sum( abs(SIGMA(:,i) - mean(SIGMA(:, i))) )); 

            
            UpperboundGaussians{iodor, iconc}{i} = gaussian1D(x, mean(MU(:, i)) + noise_mu, mean(SIGMA(:, i)) + noise_sigma);
            LowerboundGaussians{iodor, iconc}{i} = gaussian1D(x, mean(MU(:, i)) - noise_mu, mean(SIGMA(:, i)) - noise_sigma);

        end
        % Mean activation figure
        %figure(figmean)
        %hold on
        %for i = 1 : 3
        %    y1 = gaussian1D(x, mean(MU(:,i)), mean(SIGMA(:, i)));
        %    plot(x, y1, '-','Color', c(iconc,:) , 'Linewidth', 2);
        %end
        %title([num2str(odorlabels{iodor}), '_' , num2str(trials_numbers.(odorlabels{iodor}).concs(iconc))])        
    end    
end
end