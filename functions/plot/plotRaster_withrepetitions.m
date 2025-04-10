function plotRaster_withrepetitions( Rast, mycolors)
   
    Nlines = 1500; %length(Rast)*length(Rast{1}); 
    plot([0, 0], [-Nlines, 0], 'k--')
    hold on
    iiu = 0; 
    for iu = 1 : length(Rast)
        iiu = iiu + 1; 
        for irep = 1 : length(Rast{iu})
            plot(Rast{iu}{irep}, (-irep*100-iiu)*ones(size(Rast{iu}{irep})), "|", 'Color', mycolors(irep,:), 'LineWidth',0.1)
        end
    end
end

