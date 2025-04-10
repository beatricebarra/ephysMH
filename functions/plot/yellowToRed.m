function cmap = yellowToRed(n)
    cmap = [ones(n,1), linspace(1, 0, n)', zeros(n,1)]; % Red stays 1, Green decreases, Blue stays 0
    colormap(cmap);
    colorbar; % Optional: show the colorbar
end