function cmap = blackToBlue(n)
    cmap = [ zeros(n, 1), zeros(n, 1), linspace(0, 1, n)']; % Blue channel increases
    colormap(cmap);
    colorbar; % Optional: show the colorbar
end