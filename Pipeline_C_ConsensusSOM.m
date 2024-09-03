
for choose_index = 1 : size(PipelineB_Results.UseDims,1)

    % Parameters (defined in section B)
    som_dim    = PipelineB_Results.UseDims(choose_index,:) ;
    som_epochs = PipelineB_Results.UseEpochs(choose_index) ;

    clear ami ari
    % Parametri
    K_Repetitions = 100 ;
    TAU = 1 ;

    % Memorizza i risultati delle diverse esecuzioni della SOM
    results = zeros(size(reduced_data,1),K_Repetitions);

    for i = 1:K_Repetitions
        net = selforgmap(som_dim, 100, 3, UseGrid );
        net.trainParam.epochs = som_epochs;
        net.trainParam.showWindow = 0;
        net = train(net, reduced_data');
        results(:,i) = vec2ind(net(reduced_data'));
    end

    %% Consensus
    clear consensus* clusters
    consensus_matrix = zeros(size(reduced_data, 1));
    for i = 1:K_Repetitions
        for j = 1:size(reduced_data, 1)
            consensus_matrix(:,j) = consensus_matrix(:,j) + (results(:,i) == results(j, i));
        end
    end
    consensus_matrix = consensus_matrix / K_Repetitions;
    Z = linkage(squareform(1 - consensus_matrix), 'average');
    clusters = cluster(Z, 'maxclust', som_dim(1)^2 );

    % PRENDI IL CLUSTERING CHE ASSOMIGLIA DI PIU
    for i = 1:K_Repetitions
        consensusRes_ari(i) = rand_index(clusters,results(:,i), 'adjusted') ;
        consensusRes_ami(i) = ami(clusters,results(:,i)) ;
    end
    [best_ari, best_run] = max(consensusRes_ari);
    [best_ami, best_runami] = max(consensusRes_ami);
    best_clustering{choose_index} = results(:, best_run);
    som_clusters = best_clustering{choose_index} ;

    % Scrivi a schermo
    disp(' ') ;
    disp( [ 'Best Adjusted Rand Index         = ' num2str(best_ari) ' (id=' num2str(best_run) ')' ] )
    disp( [ 'Best Adjusted Mutual Information = ' num2str(best_ami) ' (id=' num2str(best_runami) ')' ] )
    disp(' ') ;

    % Visualizza ARI di tutte le runs
    figure
    histogram(consensusRes_ami)
    flnm = ['Fig_pipeC_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_HistogramAMI'] ;
    export_fig(flnm, '-jpg' , '-tif', '-r500')

    %% PERCENTUALI
    %     The code here aims to analyze the within-subject variability in SOM cluster assignments for multiple resting-state measurements.
    % The code assesses the within-subject variability in the SOM cluster assignments. It shows how many different clusters each subject's resting-state measurements are assigned to.
    % This analysis is crucial for understanding the concept of degeneracy in brain-behavior coding. Degeneracy refers to the ability of different neural configurations (clusters) to produce the same behavioral outcome.

    % By counting the unique clusters per subject, the code provides insight into the stability and consistency of the clustering results. If a subject's measurements are assigned to many different clusters, it might indicate lower stability or greater within-subject variability.
    % The histogram visualizes the distribution of cluster variability across subjects, which can help identify patterns. For example, if most subjects have measurements in only one or two clusters, it might suggest consistent clustering. Conversely, a wide distribution would indicate more variability.
    % The code snippet is important because it helps to quantify and visualize the within-subject variability in SOM cluster assignments. This information is critical for understanding the stability of the clustering results and the concept of degeneracy in brain-behavior coding. By analyzing how many different clusters each subject's resting-state measurements fall into, researchers can gain insights into the robustness of the clustering and the extent to which different neural configurations can produce similar behavioral outcomes.
    subj_percentages = zeros(numSubjects, 1);
    for i = 1 : numSubjects
        exSub1Ind = i + (numSubjects .* ((1:NumRest) - 1)); % Identify indices for the subject's resting-state measurement
        subj_percentages(i) = numel(unique(som_clusters(exSub1Ind)));   % Count the number of unique SOM clusters assigned to this subject's measurements
    end

    %     Simulate data;
    Pipeline_Z_Control_AssignationsHistogram % OUTPUT: simDataTable.NumeroClusterUnici

    % Create a figure for the histogram
    figure('Position', [100, 100, 800, 600]);
    subplot(2,1,2)
    hold on;
    histogram(simDataTable.NumeroClusterUnici, 'BinMethod', 'integer', 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'k' ,'FaceAlpha',0.6);
    ylabel('Frequency', 'FontSize', 14);
    title('Monte-Carlo Simulation', 'FontSize', 16);
    grid on; box off;
    % Add mean and standard deviation lines
    mean_val = mean(simDataTable.NumeroClusterUnici);
    xline(mean_val, '--', 'LineWidth', 2 ,'Color',[.2 .2 .2]);
    % Add a legend
    legend({'Unique Clusters', 'Mean'}, 'Location', 'eastoutside', 'FontSize', 18);
    set(gca,'FontSize',17 , 'XTick',1:NumRest,'XTickLabel',1:4)
    xlim([0 NumRest+1])

    subplot(2,1,1)
    hold on;
    histogram(subj_percentages, 'BinMethod', 'integer', 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'w' );
    ylabel('Frequency', 'FontSize', 14);
    title('Unique SOM Clusters per Subject', 'FontSize', 16);
    grid on; box off;
    % Add mean and standard deviation lines
    mean_val = mean(subj_percentages);
    xline(mean_val, '--r', 'LineWidth', 2 );
    % Add a legend
    legend({'Unique Clusters', 'Mean'}, 'Location', 'eastoutside', 'FontSize', 18);
    set(gca,'FontSize',17 , 'XTick',1:NumRest,'XTickLabel',1:4)
    xlim([0 NumRest+1])

    flnm = ['Fig_pipeC_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_subjPercentages'] ;
    export_fig(flnm, '-jpg' , '-tif', '-r500')

    %%
    % SALVA!
    save(['PipeC_Id' num2str(IdProcedure) '_Out_Consensus_Ncycles' num2str(K_Repetitions) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '.mat'] , 'best_clustering' , 'subj_percentages' )

    %% Visualizzazione delle SOM con matrici di connettività medie e informazioni aggiuntive
    %% 1) Visualizzazione NODI

    % Generare una nuova ordinazione delle regioni basata sui network
    [~, order] = sort(nets.AichaBnp);
    nets.sorted_AichaBnp = nets.AichaBnp(order); % salva indici utili
    nets.XYZ = XYZ(order,:) ;               % salva indici utili
    nets.XYZ_lab = XYZ_lab(order) ;         % salva indici utili

    figure('Position', [50, 50, 1200, 800]);

    for i = 1:som_dim(1)
        for j = 1:som_dim(2)
            neuron_index = (i-1) * som_dim(2) + j;
            cluster_indices = find(som_clusters == neuron_index);

            subplot(som_dim(1), som_dim(2), neuron_index);

            if ~isempty(cluster_indices)
                cluster_matrices = connectivity_matrices(:, :, cluster_indices);
                mean_matrix = mean(cluster_matrices, 3);

                % Riordinare la matrice di connettività media
                mean_matrix_sorted = resort_matrix(mean_matrix, order);

                % Visualizzazione della matrice di connettività media
                imagesc(mean_matrix_sorted);
                colorbar;
                axis square;

                % Informazioni aggiuntive
                num_samples = length(cluster_indices);
                clim([-.8 2]);
                title(sprintf('(X%d,Y%d) - Archetype %d (N=%d)', j, i, neuron_index, num_samples ));
                if USINGREALDATA == 1
                    colormap jet
                else
                    colormap summer
                end

            else
                title(sprintf('Archetype(%d,%d)', j, i));
                xlabel('Nessun dato');
                axis square;
            end
            set(gca,'XTick',[],'YTick',[])
        end
    end
    sgtitle('Self-Organizing Maps (SOM): Average SOM Connectomes');
    flnm = ['Fig_pipeC_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_Clusters'] ;
    export_fig(flnm, '-jpg' , '-tif', '-r500')

    Pipeline_C_ConsensusSOM_sub_plotXY

    %% 2) Visualizzazione NETWORKS
    figure('Position', [50, 50, 1200, 800]);
    usenanstd = 0 ;

    for i = 1:som_dim(1)
        for j = 1:som_dim(2)
            neuron_index = (i-1) * som_dim(2) + j;
            cluster_indices = find(som_clusters == neuron_index);

            subplot(som_dim(1), som_dim(2), neuron_index);

            if ~isempty(cluster_indices)
                cluster_matrices = connectivity_matrices(:, :, cluster_indices);
                mean_matrix = mean(cluster_matrices, 3);

                % Riordinare la matrice di connettività media
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

                % Visualizzazione della Connettività Media tra Network
                imagesc(network_connectivity);
                if usenanstd==1
                    clim([.05 .35]);
                else
                    clim([-.6 1.4]);
                end
                colorbar;
                axis square;
                num_samples = length(cluster_indices);
                title(sprintf('(X%d,Y%d) - Archetype %d (N=%d)', j, i, neuron_index, num_samples ));
                colormap jet;

                % Aggiungere le etichette degli assi
                set(gca, 'XTick', 1:num_networks, 'XTickLabel', nets.AichaBnp_labs, ...
                    'YTick', 1:num_networks, 'YTickLabel', nets.AichaBnp_labs);
                xtickangle(45);

            else
                title(sprintf('Archetype(%d,%d)', j, i));
                xlabel('Nessun dato');
                axis square;
            end
% % %             set(gca,'XTick',[],'YTick',[])
        end
    end
    sgtitle('Self-Organizing Maps (SOM): Average SOM Connectomes');

    flnm = ['Fig_pipeC_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_Clustersnets'] ;
    export_fig(flnm, '-jpg' , '-tif', '-r500')

    %% DEEPENING
    Pipeline_Cbis_ConsensusSOM_Infos   

end

%% FUNCTIONS
% Riordinare la matrice di connettività in base alla nuova ordinazione
function sorted_matrix = resort_matrix(matrix, order)
sorted_matrix = matrix(order, order);
end
