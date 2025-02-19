function [clusterCells, IDX] = cluster_units(PCMAT, clustercols)
    % PCA
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(PCMAT); 
    figure(), imagesc(COEFF)
    % KMEANS
    sil = zeros(1, 6);  
    DATA = PCMAT; 
    for K = 1 : 6
        [IDX, C, SUMD] = kmeans(DATA, K); 
        s = silhouette(DATA, IDX, 'euclidean');
        sil(K)= mean(s); 
    end
    [U, K] = max(sil); 
    [IDX, C, SUMD] = kmeans(DATA, K);
    for i = unique(IDX)'
        idx = find(IDX == i); 
        clusterCells{i} = PCMAT(idx, :); 
    end
    
    nrows = K; 
    ncols = 3; 
    cn = [1:1:K] ; 
    figure
    icol = 1; 
    selplots = icol:ncols:(nrows * ncols); 
    subplot(nrows, ncols, selplots)
    plot3(SCORE(:, 1), SCORE(:, 2), SCORE(:, 3), 'ok')
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    
    title('PCA')
    icol = 2; 
    selplots = icol:ncols:(nrows * ncols); 
    subplot(nrows, ncols, selplots)
    for i = unique(IDX)'
        idx = find(IDX == i); 
        plot(SCORE(idx, 1), SCORE(idx, 2), 'ok', 'MarkerFaceColor', clustercols{i});  
        hold on
    end
    title('Cluster Assignment (kmeans)')
    for iK = 1 : K
        subplot(nrows, ncols, ncols*iK)
        imagesc(clusterCells{iK})
        title(['Cluster ', num2str(iK)])
    end
end