function plot_session_SniffIntervals(fig, Odors, Concentration, ISI)

figure(fig)
io = 0; 
cols = ['r', 'b', 'g', 'y']; 
for o = unique(Odors)'
    
    conc_this_odor = unique(Concentration(find(Odors == o))); 
    io = io+1; 
    subplot(1, 4, io)
    
    x = []; 
    y = []; 
    allconc = unique(conc_this_odor); 
    ic = 1; 
    for c = allconc(2:end)'
        ic = ic +1; 
        x = [x; log10(c)*ones(1, length(ISI(io, ic, :)))'];
        y = [y; squeeze(ISI(io, ic, :))]; 
        scatter(log10(c)*ones(1, length(ISI(io, ic, :))), squeeze(ISI(io, ic, :)), 'ok', 'MarkerFaceColor', cols(io), 'MarkerEdgeColor', cols(io), 'MarkerFaceAlpha',0.2)
        hold on
        plot(log10(c), mean(squeeze(ISI(io, ic, :))), 'ok', 'MarkerFaceColor', cols(io), 'MarkerSize', 12)       
    end

    f = fit(x,y,'logistic', 'Lower',[0 0 -4], ...
                  'Upper',[Inf 10 -1]);
    plot(f,x,y, 'k')

    ylim([0, 2.5])
end

end