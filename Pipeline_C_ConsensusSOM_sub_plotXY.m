
%% THIS WORKS ONLY FO 2x2 MATRICES!

%% 1) X1
figure('Position', [50, 50, 700, 450]);

% X1: A1 e A3
neuron_index1 = 1 ; % A1
neuron_index2 = 3 ; % A3
% estrai indici e calcola matrice
cluster_indices = [ find(som_clusters == neuron_index1) ; find(som_clusters == neuron_index2) ];
cluster_matrices = connectivity_matrices(:, :, cluster_indices);
mean_matrix = mean(cluster_matrices, 3);
% Riordinare la matrice di connettività media
mean_matrix_sorted = resort_matrix(mean_matrix, order);

% Visualizzazione della matrice di connettività media
imagesc(mean_matrix_sorted);
colorbar;
axis square;
% Informazioni aggiuntive
clim([-.8 2]);
title('X1: A1 & A3');
colormap jet
set(gca,'XTick',[],'YTick',[])

% Esporta
flnm = ['Fig_pipeC_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_Clusters_X1'] ;
export_fig(flnm, '-jpg' , '-tif', '-r500')

%% 2) X2
figure('Position', [50, 50, 700, 450]);

% X2: A2 e A4
neuron_index1 = 2 ; % A2
neuron_index2 = 4 ; % A4
% estrai indici e calcola matrice
cluster_indices = [ find(som_clusters == neuron_index1) ; find(som_clusters == neuron_index2) ];
cluster_matrices = connectivity_matrices(:, :, cluster_indices);
mean_matrix = mean(cluster_matrices, 3);
% Riordinare la matrice di connettività media
mean_matrix_sorted = resort_matrix(mean_matrix, order);

% Visualizzazione della matrice di connettività media
imagesc(mean_matrix_sorted);
colorbar;
axis square;
% Informazioni aggiuntive
clim([-.8 2]);
title('X2: A2 & A4');
colormap jet
set(gca,'XTick',[],'YTick',[])

% Esporta
flnm = ['Fig_pipeC_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_Clusters_X2'] ;
export_fig(flnm, '-jpg' , '-tif', '-r500')

%% 3) Y1
figure('Position', [50, 50, 700, 450]);

% Y1: A1 e A2
neuron_index1 = 1 ; % A1
neuron_index2 = 2 ; % A2
% estrai indici e calcola matrice
cluster_indices = [ find(som_clusters == neuron_index1) ; find(som_clusters == neuron_index2) ];
cluster_matrices = connectivity_matrices(:, :, cluster_indices);
mean_matrix = mean(cluster_matrices, 3);
% Riordinare la matrice di connettività media
mean_matrix_sorted = resort_matrix(mean_matrix, order);

% Visualizzazione della matrice di connettività media
imagesc(mean_matrix_sorted);
colorbar;
axis square;
% Informazioni aggiuntive
clim([-.8 2]);
title('Y1: A1 & A2');
colormap jet
set(gca,'XTick',[],'YTick',[])

% Esporta
flnm = ['Fig_pipeC_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_Clusters_Y1'] ;
export_fig(flnm, '-jpg' , '-tif', '-r500')

%% 4) Y2
figure('Position', [50, 50, 700, 450]);

% Y2: A3 e A4
neuron_index1 = 3 ; % A3
neuron_index2 = 4 ; % A4
% estrai indici e calcola matrice
cluster_indices = [ find(som_clusters == neuron_index1) ; find(som_clusters == neuron_index2) ];
cluster_matrices = connectivity_matrices(:, :, cluster_indices);
mean_matrix = mean(cluster_matrices, 3);
% Riordinare la matrice di connettività media
mean_matrix_sorted = resort_matrix(mean_matrix, order);

% Visualizzazione della matrice di connettività media
imagesc(mean_matrix_sorted);
colorbar;
axis square;
% Informazioni aggiuntive
clim([-.8 2]);
title('Y2: A3 & A4');
colormap jet
set(gca,'XTick',[],'YTick',[])

% Esporta
flnm = ['Fig_pipeC_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_Clusters_Y2'] ;
export_fig(flnm, '-jpg' , '-tif', '-r500')

%%
% Riordinare la matrice di connettività in base alla nuova ordinazione
function sorted_matrix = resort_matrix(matrix, order)
sorted_matrix = matrix(order, order);
end
