function [ISI] = plot_session_respiration(Odors, Concentration, FVO, FVC, fSniff, SniffTime, Prex)
    
figure
    io = 0; 
    for o = unique(Odors)'
        conc_this_odor = unique(Concentration(find(Odors == o))); 
        ic = 1;
        io = io+1; 
        for c = unique(conc_this_odor)'
            
            subplot(length(unique(Odors)), length(unique(conc_this_odor)), double((o-5)*length(unique(conc_this_odor))+ic))
            idx_trials = find(Odors == o & Concentration == c); 
            ii = 0;
            iii = 0; 
            vertical_offset = 0.3; 
            for itrial =  idx_trials'
                iii = iii +1; 
                ii = ii -vertical_offset;
                start_time = FVO(itrial) - 10;
                end_time = FVO(itrial) + 6;
                start_idx = searchsorted(SniffTime, start_time);
                end_idx = searchsorted(SniffTime, end_time);
                % Extract respiration traces
                respiration_trace_window = fSniff(start_idx:end_idx)+ii; 
                postfvoinh = find(Prex>FVO(itrial)); 
                inh1 = Prex(postfvoinh(1)); 
                inh2 = Prex(postfvoinh(2)); 
                % Save duration of inter-sniff-duration
                ISI(io, ic, iii) = inh2 - inh1;
                % Normalize time axis relative to FVO
                trace_times = SniffTime(start_idx:end_idx) - FVO(itrial); 
                
                hold on
                fill([0, FVC(itrial)-FVO(itrial), FVC(itrial)-FVO(itrial), 0 ],...
                    [mean(respiration_trace_window)-vertical_offset/2, mean(respiration_trace_window)-vertical_offset/2, ...
                    mean(respiration_trace_window)+vertical_offset/2, mean(respiration_trace_window)+vertical_offset/2],  'b', 'FaceAlpha', 0.4, 'EdgeColor', 'b', 'EdgeAlpha', 0.2)
            
                fill([inh1-FVO(itrial), inh2-FVO(itrial), inh2-FVO(itrial), inh1-FVO(itrial) ],...
                    [mean(respiration_trace_window)-vertical_offset/2, mean(respiration_trace_window)-vertical_offset/2, ...
                    mean(respiration_trace_window)+vertical_offset/2, mean(respiration_trace_window)+vertical_offset/2],  [0.8, 0.8, 0.8], 'FaceAlpha', 1, 'EdgeColor', [0.8, 0.8, 0.8], 'EdgeAlpha', 0.5)
                
                plot(trace_times, respiration_trace_window,'k', 'Linewidth', 1.5)
                title(['Vial = ', num2str(o), ',  Conc = ', num2str(c)])
                text(trace_times(end) + 0.2, mean(respiration_trace_window), ['Trial # ', num2str(itrial)])
                percent_complete = (iii / length(idx_trials)) * 100;
                disp([' ----   ', num2str(percent_complete), ' % of Sniff plotting run   ----'])
            end
    
            ic = ic +1; 
            xlim([-0.5,2])
        end
    end

end