function cmap = blackToRed(n)
    cmap = [linspace(0, 1, n)', zeros(n, 1), zeros(n, 1)]; % Red channel increases
    colormap(cmap);
    colorbar; % Optional: show the colorbar
end