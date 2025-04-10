function cmap = blue_white_red(n)
    if nargin < 1
        n = 256; % Default colormap size
    end
    
    % Define the key colors: blue, white, red
    c1 = [0, 0, 1];   % Blue
    c2 = [1, 1, 1];   % White
    c3 = [1, 0, 0];   % Red
    
    % Interpolation points
    x = linspace(0, 1, n);
    xq = [0, 0.5, 1];
    
    % Interpolate each channel
    r = interp1(xq, [c1(1), c2(1), c3(1)], x);
    g = interp1(xq, [c1(2), c2(2), c3(2)], x);
    b = interp1(xq, [c1(3), c2(3), c3(3)], x);
    
    % Combine into colormap
    cmap = [r(:), g(:), b(:)];
end
