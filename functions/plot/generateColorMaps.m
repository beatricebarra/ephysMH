function colormaps = generateColorMaps(N)
% generateColorMaps Generates 4 colormaps: red, blue, yellow, and green
% Each colormap contains N colors from light to dark.
% Returns a struct with fields: red, blue, yellow, green

    if nargin < 1
        N = 100; % default number of colors
    end

    % Light to dark mapping: linearly interpolate from light color to pure color
    colormaps.red    = interpolateColor([1, 0.8, 0.8], [0.5, 0, 0], N);
    colormaps.blue   = interpolateColor([0.8, 0.8, 1], [0, 0, 0.5], N);
    colormaps.yellow = interpolateColor([1, 1, 0.8], [0.5, 0.5, 0], N);
    colormaps.green  = interpolateColor([0.8, 1, 0.8], [0, 0.5, 0], N);

end

function cmap = interpolateColor(lightRGB, darkRGB, N)
% Helper function to interpolate between light and dark RGB colors
    cmap = zeros(N, 3);
    for i = 1:3
        cmap(:, i) = linspace(lightRGB(i), darkRGB(i), N);
    end
end