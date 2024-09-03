
numClusters = som_dim(1)^2 ;

% Inizializzare una matrice per memorizzare le differenze medie per ogni neurone
differences = zeros(Nregions, Nregions, numClusters);

% Calcolare la matrice di connettività media per ciascun neurone (cluster)
mean_matrices = zeros(Nregions, Nregions, numClusters);
for k = 1:numClusters
    cluster_indices = find(som_clusters == k);
    if ~isempty(cluster_indices)
        cluster_matrices = connectivity_matrices(:, :, cluster_indices);
        mean_matrices(:, :, k) = mean(cluster_matrices, 3);
    end
end

% Calcolare le differenze tra ogni neurone e tutti gli altri
for k = 1:numClusters
    for j = 1:numClusters
        if k ~= j
            differences(:, :, k) = differences(:, :, k) + abs(mean_matrices(:, :, k) - mean_matrices(:, :, j));
        end
    end
end

% Media delle differenze per ogni connessione
mean_differences = differences / (numClusters - 1);

% Identificare le connessioni che cambiano di più per ogni neurone
significant_changes = false(Nregions, Nregions, numClusters);
significant_stables = NaN(Nregions, Nregions, numClusters);
for k = 1:numClusters
    tempMeanDiff = mean_differences(:, :, k) ;
    tempMeanDiffLinear = tempMeanDiff(lowerTri>0) ;
    threshold = prctile(tempMeanDiffLinear(:), 95); % Consideriamo il 95° percentile delle differenze per ogni neurone
    threshold2 = prctile(tempMeanDiffLinear(:), 5); % Consideriamo il 95° percentile delle differenze per ogni neurone
    significant_changes(:, :, k) = tempMeanDiff > threshold;
    significant_stables(:, :, k) = tempMeanDiff < threshold2 & tempMeanDiff > 0;
end

% Impostare i limiti di colore per la visualizzazione
clim_min = min(mean_differences(:));
clim_max = max(mean_differences(:));

%% Visualizzazione delle differenze medie (NON UTILE)
figure('Position', [50, 50, 1200, 800]);
for k = 1:numClusters
    subplot(som_dim(1), som_dim(2), k);
    mean_matrix = mean_differences(:, :, k) ;
    mean_matrix_sorted = resort_matrix(mean_matrix, order);
    imagesc(mean_matrix_sorted);
    colorbar;
    title(sprintf('Average Differences - Archetype %d', k));
    axis square;
    clim([clim_min clim_max]);
    colormap jet;
end
sgtitle('Average Archetypal Differences (POCO UTILE)');
flnm = ['Fig_pipeCbis_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_AvgClustDiffs'] ;
export_fig(flnm, '-jpg' , '-tif', '-r500') 

% Visualizzazione delle connessioni che cambiano di più
figure('Position', [50, 50, 1200, 800]);
for k = 1:numClusters
    subplot(som_dim(1), som_dim(2), k);
    mean_matrix = significant_changes(:, :, k) ;
    mean_matrix_sorted = resort_matrix(mean_matrix, order);
    imagesc(mean_matrix_sorted);
    colorbar;
    title(sprintf('Modulated Connections - Archetype %d', k));
    axis square;
    colormap summer;
    set(gca,'XTick',[],'YTick',[])
end
sgtitle('5% most modulated connections across Archetypes');
flnm = ['Fig_pipeCbis_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_ConnModulated'] ;
export_fig(flnm, '-jpg' , '-tif', '-r500') 

% Visualizzazione delle connessioni che cambiano di meno
figure('Position', [50, 50, 1200, 800]);
for k = 1:numClusters
    subplot(som_dim(1), som_dim(2), k);
    mean_matrix = significant_stables(:, :, k);
    mean_matrix_sorted = resort_matrix(mean_matrix, order);
    imagesc(mean_matrix_sorted);
    colorbar;
    title(sprintf('Stable Connections - Archetype %d', k));
    axis square;
    colormap summer;
    set(gca,'XTick',[],'YTick',[])
end
sgtitle('5% most stable connections across Archetypes');
flnm = ['Fig_pipeCbis_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_ConnStable'] ;
export_fig(flnm, '-jpg' , '-tif', '-r500') 

%%
% Analisi delle differenze tra coppie di neuroni
pairwise_differences = zeros(Nregions, Nregions, numClusters, numClusters);

for k = 1:numClusters
    for j = 1:numClusters
        if k ~= j
            pairwise_differences(:, :, k, j) = (mean_matrices(:, :, k) - mean_matrices(:, :, j));
        end
    end
end

% Visualizzazione delle differenze tra coppie di neuroni
figure('Position', [50, 50, 1200, 800]);
for k = 1:numClusters
    for j = 1:numClusters
        if k ~= j
            subplot(numClusters, numClusters, (k-1)*numClusters + j);
            mean_matrix = pairwise_differences(:, :, k, j) ;
            mean_matrix_sorted = resort_matrix(mean_matrix, order);
            imagesc(mean_matrix_sorted);
            colorbar;
            title(sprintf('A%d - A%d', k, j));
            axis square;
            clim([clim_min clim_max]);
            set(gca,'XTick',[],'Ytick',[])
            colormap jet;
        end
    end
end
sgtitle('Differences across couples of archetypes');
flnm = ['Fig_pipeCbis_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_ClustDiffs'] ;
export_fig(flnm, '-jpg' , '-tif', '-r500') 

%% NOW NETWORKS

%% Visualizzazione delle differenze medie (NON UTILE)
figure('Position', [50, 50, 1200, 800]);
for k = 1:numClusters
    subplot(som_dim(1), som_dim(2), k);

    mean_matrix = mean_differences(:, :, k) ;
    mean_matrix_sorted = resort_matrix(mean_matrix, order);

    % Calcolo FC media
    unique_networks = unique(nets.AichaBnp);
    num_networks = length(unique_networks);
    network_connectivity = zeros(num_networks, num_networks);

    for m = 1:num_networks
        for n = 1:num_networks
            regions_m = find(nets.sorted_AichaBnp == unique_networks(m));
            regions_n = find(nets.sorted_AichaBnp == unique_networks(n));
            conn_values = mean_matrix_sorted(regions_m, regions_n);
            if usenanstd==1
                network_connectivity(m, n) = nanstd(conn_values(:));
            else
                network_connectivity(m, n) = nanmean(conn_values(:));
            end
        end
    end

    imagesc(network_connectivity);
    colorbar;
    title(sprintf('Average Differences - Archetype %d', k));
    axis square;
    clim([clim_min clim_max]);
    clim([0 0.5])
    colormap jet;
end
sgtitle('Average Archetypal Differences (POCO UTILE)');
flnm = ['Fig_pipeCbis_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_NET_AvgClustDiffs'] ;
export_fig(flnm, '-jpg' , '-tif', '-r500') 

% Visualizzazione delle connessioni che cambiano di più
figure('Position', [50, 50, 1200, 800]);
for k = 1:numClusters
    subplot(som_dim(1), som_dim(2), k);

    mean_matrix = significant_changes(:, :, k) ;
    mean_matrix_sorted = resort_matrix(mean_matrix, order);
    
    % Calcolo FC media
    unique_networks = unique(nets.AichaBnp);
    num_networks = length(unique_networks);
    network_connectivity = zeros(num_networks, num_networks);

    for m = 1:num_networks
        for n = 1:num_networks
            regions_m = find(nets.sorted_AichaBnp == unique_networks(m));
            regions_n = find(nets.sorted_AichaBnp == unique_networks(n));
            conn_values = mean_matrix_sorted(regions_m, regions_n);
            if usenanstd==1
                network_connectivity(m, n) = nanstd(conn_values(:));
            else
                network_connectivity(m, n) = nanmean(conn_values(:));
            end
        end
    end
    
    imagesc(network_connectivity);
    colorbar;
    title(sprintf('Modulated Connections - Archetype %d', k));
    axis square;
    colormap summer;
    set(gca,'XTick',1:numel(nets.AichaBnp_labs),'YTick',1:numel(nets.AichaBnp_labs),'XTickLabel',nets.AichaBnp_labs,'YTickLabel',nets.AichaBnp_labs,'XTickLabelRotation',45)
end
sgtitle('5% most modulated connections across Archetypes');
flnm = ['Fig_pipeCbis_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_NET_ConnModulated'] ;
export_fig(flnm, '-jpg' , '-tif', '-r500') 

% Visualizzazione delle connessioni che cambiano di meno
figure('Position', [50, 50, 1200, 800]);
for k = 1:numClusters
    subplot(som_dim(1), som_dim(2), k);

    mean_matrix = significant_stables(:, :, k);
    mean_matrix_sorted = resort_matrix(mean_matrix, order);
    
    % Calcolo FC media
    unique_networks = unique(nets.AichaBnp);
    num_networks = length(unique_networks);
    network_connectivity = zeros(num_networks, num_networks);

    for m = 1:num_networks
        for n = 1:num_networks
            regions_m = find(nets.sorted_AichaBnp == unique_networks(m));
            regions_n = find(nets.sorted_AichaBnp == unique_networks(n));
            conn_values = mean_matrix_sorted(regions_m, regions_n);
            if usenanstd==1
                network_connectivity(m, n) = nanstd(conn_values(:));
            else
                network_connectivity(m, n) = nanmean(conn_values(:));
            end
        end
    end
    
    imagesc(network_connectivity);
    colorbar;
    title(sprintf('Stable Connections - Archetype %d', k));
    axis square;
    colormap summer;
    set(gca,'XTick',[],'YTick',[])
end
sgtitle('5% most stable connections across Archetypes');
flnm = ['Fig_pipeCbis_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_NET_ConnStable'] ;
export_fig(flnm, '-jpg' , '-tif', '-r500') 

%% FUNCTIONS

% Funzione per mappare i valori del triangolo inferiore in una matrice completa
function full_matrix = map_to_full_matrix(vector, matrix_size)
    full_matrix = zeros(matrix_size);
    index = 1;
    for i = 1:matrix_size
        for j = 1:i-1
            full_matrix(i, j) = vector(index);
            full_matrix(j, i) = vector(index); % Assicurare la simmetria
            index = index + 1;
        end
    end
end

%% FUNCTIONS
% Funzione per Riordinare la matrice di connettività in base alla nuova ordinazione
function sorted_matrix = resort_matrix(matrix, order)
    sorted_matrix = matrix(order, order);
end
