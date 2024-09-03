
clear
close all hidden
clc

IdProcedure = 99 ;

% load best clustering
load('PipeC_Id99_Out_Consensus_Ncycles100_Dim2_epochs200.mat')
load('C:\Users\Hp\Desktop\PostDoc2023\ResearchTopic_BRAINBEHAVIOR_Frontiers\BrainGame_Reloaded\BrainGame_2024\RealData_MultiRest\RESULTS_DEFINITIVI\PipeB_Id99_Out_Parametrization_Ncycles4.mat')
load('../XYZ_AI.mat')

% other
addpath(genpath('C:\Users\Hp\Documents\MATLAB\BrainNetViewer\'))
addpath(genpath('C:\Users\Hp\Documents\MATLAB\BCT\'))
addpath(genpath('C:\Users\Hp\Documents\MATLAB\Export_Fig\'))
addpath(genpath('C:\Users\Hp\Documents\MATLAB\othercolor\'))
addpath(genpath('C:\Users\Hp\Documents\MATLAB\subplot_tight\'))

%% PRIORS
thr_absconnDiffs = 0.1 ; % soglia assoluta di differenze FC
thr = 99.9 ;              % percentuale di FC differenti ritenute (dopo applicazione soglia linea precedente)

choose_index = 1 ;
som_clusters = best_clustering{choose_index} ;
som_dim    = [2 2] ;
som_epochs = 100 ;
numClusters = som_dim(1)^2 ;
NumRest = 4 ;
Nregions = 344 ;
numSubjects = 496 ;        % Choose number of subjects

%% LOAD FC
load('FC_ALL_496subjgood_4rest.mat')                % Import Brain

connectivity_matrices = [] ;
% Prendo solo la corteccia (1:344)
for r = 1 : NumRest
    preconnectivity_matrices = FC.z(1:Nregions,1:Nregions,1:numSubjects,r);
    connectivity_matrices = cat ( 3 ,connectivity_matrices , preconnectivity_matrices) ;
end
connectivity_matrices (abs(connectivity_matrices)>4) = NaN ;
for n = 1:numSubjects*NumRest
    lowerTri = tril(connectivity_matrices(:,:,n), -1);
    flattened_matrices(n, :) = lowerTri(lowerTri ~= 0);
end

load('ConersionMat_NoFC.mat'); % Uso ConversionMat
AichaBnp = zeros (Nregions,1) ;
% per ogni riga (regione)
for i = 1 : Nregions
    % attribuisci una sola colonna (network)
    [~,id] = max(ConversionMat(i,:)) ;
    AichaBnp(i) = id ;
end

% Load and rearrange Labels
load("cfrnets_7_16.mat") % Labels are in Ganv.labs16
TheNetOrder = [ 6 , 1 , 4 , 3 , 14 , 13 , ...
    5 , 10 , 12 , 16 , 15 , ...
    2 , 7 , 11 , 8 , 9 ] ;
nets.PreviousOrder = fliplr(TheNetOrder) ;

% Reassign XYZ
dep = AichaBnp ;
for n = 1 : numel(nets.PreviousOrder)
    dep(AichaBnp==nets.PreviousOrder(n)) = n ;
end
nets.AichaBnp = dep ;
nets.AichaBnp_labs = Ganv.labs16(nets.PreviousOrder);
[~, order] = sort(nets.AichaBnp);
nets.sorted_AichaBnp = nets.AichaBnp(order); % salva indici utili
nets.XYZ = XYZ(order,:) ;               % salva indici utili
nets.XYZ_lab = XYZ_lab(order) ;         % salva indici utili

%% FIGURA CLUSTDIFFS

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
% NOTA CHE QUESTE SONO GIA CALCOLATE IN VALORE ASSOLUTO!
clear ttest*
for k = 1:numClusters
    for j = 1:numClusters
        if k ~= j
            differences(:, :, k) = differences(:, :, k) + abs(mean_matrices(:, :, k) - mean_matrices(:, :, j));

%             LA SEGUENTE PARTE NON LA FACCIO PIU PERCHE VIENE TUTTO
%             SIGNIFICATIVO, non ha senso
% % % % %             % CALCOLA DIFFERENZE SIGNIFICATIVE
% % % % %             Kcluster_indices = find(som_clusters == k);
% % % % %             Kcluster_matrices = connectivity_matrices(:, :, Kcluster_indices);
% % % % %             Jcluster_indices = find(som_clusters == j);
% % % % %             Jcluster_matrices = connectivity_matrices(:, :, Jcluster_indices);
% % % % %             %             per ogni connessione:
% % % % %             for region1 = 1 : Nregions
% % % % %                 for region2 = 1 : Nregions
% % % % %                     if region1~=region2
% % % % % 
% % % % %                         fc_k = Kcluster_matrices(region1, region2, : );
% % % % %                         fc_j = Jcluster_matrices(region1, region2, : );
% % % % %                         [ttest_h,ttest_p,ttest_ci,ttest_stats] = ttest2  ( fc_k(:) , fc_j(:) ) ;
% % % % % 
% % % % %                         ttest_matrix(region1,region2) = ttest_h;
% % % % % 
% % % % %                     end % if they're different
% % % % %                 end % region2
% % % % %             end % region1

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

%% Figure

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

for k = 1 : numClusters
    for j = 1 : numClusters
        if k ~= j
            %             subplot(numClusters, numClusters, (k-1)*numClusters + j);
            figure('Position', [50, 50, 1200, 800]);
            mean_matrix = pairwise_differences(:, :, k, j) ;
            % imposto a zero tutto ciò che è negativo per mantenere  "direzionalità" nelle diagonali
            mean_matrix(mean_matrix<0) = 0 ;
            mean_matrix_sorted = resort_matrix(mean_matrix, order);
            imagesc(mean_matrix_sorted);
            colorbar;
            title(sprintf('A%d - A%d', k, j));
            axis square;
            clim([clim_min clim_max]);
            set(gca,'XTick',[],'Ytick',[])
            colormap jet;
            sgtitle('Differences across couples of archetypes');
            flnm = ['Fig_E_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_ClustDiffs_k' num2str(k) '_j' num2str(j)] ;
            export_fig(flnm, '-jpg' , '-tif', '-r500' )

            %% BRAIN NET VIEWER SECTION!!
            %% PREPARE DATA
            % Mantain only strongest modulations
            sub_mean_matrix_sorted = mean_matrix_sorted ;
            sub_mean_matrix_sorted(sub_mean_matrix_sorted<thr_absconnDiffs) = NaN ;
            sub_XYZ = nets.XYZ ;
            sub_XYZ_lab = nets.XYZ_lab ;
            selecttriangle = tril(ones(Nregions,Nregions),-1) ;
            thr_conn = prctile(sub_mean_matrix_sorted(selecttriangle==1),thr) ;
            sub_mean_matrix_sorted(sub_mean_matrix_sorted<thr_conn) = NaN ;

            % rimuovi nodi non necessari
            idUnnecessary = find ( sum(isnan(sub_mean_matrix_sorted),2) == Nregions ) ;
            sub_mean_matrix_sorted( : , idUnnecessary ) = [] ;
            sub_mean_matrix_sorted( idUnnecessary , : ) = [] ;
            sub_XYZ (idUnnecessary,:) = [] ;
            sub_XYZ_lab (idUnnecessary) = [] ;

            %% PLOT

            % Crea il file .node
            node_file = 'nodes.node';
            fid = fopen(node_file, 'w');
            for i = 1:size(sub_XYZ, 1)
                fprintf(fid, '%f %f %f %d %d %s\n', sub_XYZ(i, 1), sub_XYZ(i, 2), sub_XYZ(i, 3), 1, 1, sub_XYZ_lab{i}); % X Y Z COLOR SIZE LABEL
            end
            fclose(fid);

            % Crea il file .edge
            edge_file = 'edges.edge';
            dlmwrite(edge_file, sub_mean_matrix_sorted, 'delimiter', '\t');

            disp(['Number of regions = ' num2str(size(sub_XYZ,1))])
            disp(['Number of connections = ' num2str(sum(sum(isnan(sub_mean_matrix_sorted)==0))/2)])

            % Carica i file in Brain Net Viewer
            flnm = ['Fig_E_brain_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_ClustDiffs_k' num2str(k) '_j' num2str(j)] ;
            BrainNet_MapCfg('C:\Users\Hp\Documents\MATLAB/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152.nv', ...
                node_file, edge_file, 'Pipeline_E_Brain_BNVSettings.mat' , [flnm '.jpg'] );

            % SAVE RESULTS
            brains.node{k,j} = table (sub_XYZ , sub_XYZ, sub_XYZ, sub_XYZ_lab, ...
                'VariableNames', {'X' 'Y' 'Z' 'LABEL' }) ;
            brains.edge{k,j} = sub_mean_matrix_sorted ;

        end % if k and k are different
    end % j
end % k


%% NEW
%% Visualizzazione delle differenze medie tra livelli delle dimensioni X e Y

% Parametri SOM (da adattare in base ai tuoi dati)
som_dim = [2, 2];  % Dimensioni della SOM 2x2
Nregions = size(pairwise_differences, 1);  % Numero di regioni

% Ciclo sulle dimensioni
for dim = 1:2  % Dimensioni X e Y
    if dim == 1
        % Differenze lungo la dimensione X (X1 vs X2)
        A1a = 1; A1b = 3;   % Livelli X1 
        A2a = 2; A2b = 4;   % Livelli X2
        title_dim = 'X';
    else
        % Differenze lungo la dimensione Y (Y1 vs Y2)
        A1a = 1; A1b = 2;   % Livelli Y1 
        A2a = 3; A2b = 4;   % Livelli Y2
        title_dim = 'Y';
    end
    
    % Calcolo delle differenze medie
    mean_matrix_1 = mean ( mean(pairwise_differences(:, :, [A1a A1b], :), 4) , 3);
    mean_matrix_2 = mean ( mean(pairwise_differences(:, :, [A2a A2b], :), 4) , 3);
    
    diff_1_2 = mean_matrix_1 - mean_matrix_2;
    diff_2_1 = mean_matrix_2 - mean_matrix_1;
    
    % Visualizzazione e salvataggio delle matrici
    for i = 1:2
        if i == 1
            mean_matrix_sorted = diff_1_2;
            title_str = sprintf('%s1_%s2', title_dim, title_dim);
        else
            mean_matrix_sorted = diff_2_1;
            title_str = sprintf('%s2_%s1', title_dim, title_dim);
        end
        
        mean_matrix_sorted(mean_matrix_sorted<0) = 0;
        mean_matrix_sorted = resort_matrix(mean_matrix_sorted, order);
        
        figure('Position', [50, 50, 1200, 800]);
        imagesc(mean_matrix_sorted);
        colorbar;
        if i == 1
            title(sprintf('%s1 > %s2', title_dim, title_dim));
        else
            title(sprintf('%s2 > %s1', title_dim, title_dim));
        end
        axis square;
        clim([0 .6]);
        set(gca,'XTick',[],'Ytick',[]);
        colormap jet;
        sgtitle(['Mean Differences across ', title_dim, ' Dimension']);
        flnm = ['Fig_MeanDiff_', title_dim, '_', title_str, '_Id', num2str(IdProcedure), '_Dim', num2str(som_dim(1)), '_epochs', num2str(som_epochs)];
        export_fig(flnm, '-jpg', '-tif', '-r400');
        
        %% BRAIN NET VIEWER SECTION!!
        %% PREPARE DATA
        sub_mean_matrix_sorted = mean_matrix_sorted;
        sub_mean_matrix_sorted(sub_mean_matrix_sorted<thr_absconnDiffs) = NaN;
        sub_XYZ = nets.XYZ;
        sub_XYZ_lab = nets.XYZ_lab;
        selecttriangle = tril(ones(Nregions,Nregions),-1);
        thr_conn = prctile(sub_mean_matrix_sorted(selecttriangle==1), thr);
        sub_mean_matrix_sorted(sub_mean_matrix_sorted<thr_conn) = NaN;
        
        % Rimuovi nodi non necessari
        idUnnecessary = find(sum(isnan(sub_mean_matrix_sorted), 2) == Nregions);
        sub_mean_matrix_sorted(:, idUnnecessary) = [];
        sub_mean_matrix_sorted(idUnnecessary, :) = [];
        sub_XYZ(idUnnecessary,:) = [];
        sub_XYZ_lab(idUnnecessary) = [];
        
        %% PLOT
        % Crea il file .node
        node_file = 'nodes.node';
        fid = fopen(node_file, 'w');
        for n = 1:size(sub_XYZ, 1)
            fprintf(fid, '%f %f %f %d %d %s\n', sub_XYZ(n, 1), sub_XYZ(n, 2), sub_XYZ(n, 3), 1, 1, sub_XYZ_lab{n});
        end
        fclose(fid);
        
        % Crea il file .edge
        edge_file = 'edges.edge';
        dlmwrite(edge_file, sub_mean_matrix_sorted, 'delimiter', '\t');
        
        disp(['Number of regions = ' num2str(size(sub_XYZ, 1))]);
        disp(['Number of connections = ' num2str(sum(sum(~isnan(sub_mean_matrix_sorted))) / 2)]);
        
        % Carica i file in Brain Net Viewer
        flnm = ['Fig_MeanDiff_Brain_', title_dim, '_', title_str, '_Id', num2str(IdProcedure), '_Dim', num2str(som_dim(1)), '_epochs', num2str(som_epochs)];
        BrainNet_MapCfg('C:\Users\Hp\Documents\MATLAB/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152.nv', ...
            node_file, edge_file, 'Pipeline_E_Brain_BNVSettings.mat', [flnm, '.jpg']);
        
        % Salva i risultati
        brains.node_mean{dim, i} = table(sub_XYZ, sub_XYZ, sub_XYZ, sub_XYZ_lab, 'VariableNames', {'X', 'Y', 'Z', 'LABEL'});
        brains.edge_mean{dim, i} = sub_mean_matrix_sorted;
    end
end

%% FUNCTIONS
% Funzione per Riordinare la matrice di connettività in base alla nuova ordinazione
function sorted_matrix = resort_matrix(matrix, order)
sorted_matrix = matrix(order, order);
end
