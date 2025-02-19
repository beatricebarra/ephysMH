function [fig] = plotRasters_withrepetitions(fig, plot_geometry, Rast, odor_idx, conc_idx, mycolors)

    % plot_geometry corresponds to Nrows, Ncols, Orows, Ocols (O = offset)
    ii = 0; 
    figure(fig)
    for i = odor_idx
        for j = conc_idx
        
            ii = ii +1; 
            subplot(length(odor_idx), length(conc_idx), ii)
            Nlines = length(Rast{i,j})*length(Rast{i,j}{iu}); 
            plot([0, 0], [-Nlines, 0], 'k--')
            hold on
            iiu = 0; 
            for iu = 1 : length(Rast{i,j})
                iiu = iiu + 1; 
                for irep = 1 : length(Rast{i,j}{iu})
                    plot(Rast{i,j}{iu}{irep}, (-irep*100-iiu)*ones(size(Rast{i,j}{iu}{irep})), "|", 'Color', mycolors(irep,:))
                end
            end
        end
    end

end

