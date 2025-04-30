function realI = computeIntensity(intensity_file, odor_concentrations)
%     params = mean(xlsread(intensity_file)); 
%     
%     %pure_ppm_odors = {'ethyltiglate':5617.105263, '2-heptanone': 5065.789474, 'ethylbutyrate':16842.10526, 'acetophenone': 526.16}
%     % ET, 2-HEP, EB, ACE : this is the order of parameters in the perceived
%     % intensity curves
%     % odorlabels = {'ethyltiglate', 'ethylbutyrate', 'acetophenone', 'heptanone'}; 
%     % This is the order in the ephys data labels and colors
%     parammat = [];
%     parammat(1,:) = params(1:3); % Ethyltiglate
%     parammat(2,:) = params(7:9); % Ethyl butyrate
%     parammat(3,:) = params(10:12); % Acetophenone
%     parammat(4,:) = params(4:6); % 2-Heptanone
%     ppm_col = [5617.105263; 16842.10 ; 526.16; 5065.789474];
%     parammat = [parammat, ppm_col]; 
    % Colors
    int_curves_cols = ['r', 'b', 'y', 'g']; 
    
    
%     figure
%     hold on
%     for io = 1 : 4
%         x = logspace(-1, 6, 100); 
%         x = log10(x); 
%         y = psychocurve(x, parammat(io, 1),  parammat(io, 2), parammat(io, 3)); 
%         
%         plot(x, y, int_curves_cols(io),  'Linewidth', 2) 
%         %set(gca,'Xscale', 'log')
%     end
    for ii = 1 : size(odor_concentrations(:,2), 1)
        % Injecting noise into the system
        allparams = xlsread(intensity_file); 
        idx_random = floor(rand(1)*size(allparams, 1)); 
        if idx_random == 0
            idx_random = 1; 
        end
        %params = allparams(idx_random,:);
        params = mean(allparams); 
        parammat = [];
        parammat(1,:) = params(1:3); % Ethyltiglate
        parammat(2,:) = params(7:9); % Ethyl butyrate
        parammat(3,:) = params(10:12); % Acetophenone
        parammat(4,:) = params(4:6); % 2-Heptanone
        ppm_col = [5617.105263; 16842.10 ; 526.16; 5065.789474];
        parammat = [parammat, ppm_col]; 

        % convert concentration ot intensity
        io = odor_concentrations(ii, 1); 
        realI(ii) = psychocurve(log10(odor_concentrations(ii, 2)*parammat(io, 4)), parammat(io, 1), parammat(io, 2), parammat(io, 3)); 
        %plot(log10(odor_concentrations(ii, 2)*parammat(io, 4)), realI(ii), 'o', 'MarkerFaceColor', int_curves_cols(io)) 
    
    end
end