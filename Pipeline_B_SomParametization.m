
if exist("IdProcedure")==0
    IdProcedure = 99 ;
end

%% INITIALIZE
clear ami* rand* adjust* mean_quantization_error legend_entries scDim
% clc

% Cell arrays per salvare i risultati
results_ami = cell(size(som_dims, 1), length(epoch_values));
results_ari = cell(size(som_dims, 1), length(epoch_values));

%% CORE PROCESS

for dim_idx = 1:size(som_dims, 1)
    for epoch_idx = 1:length(epoch_values)
        som_dim = som_dims(dim_idx, :);
        epochs = epoch_values(epoch_idx);
        som_clusters_all = zeros(num_runs, size(reduced_data, 1));

        % Eseguire il clustering SOM per num_runs volte
        for run = 1:num_runs
            % clustering
            net = selforgmap(som_dim, 100, 3, UseGrid );
            net.trainParam.epochs = epochs;
            net.trainParam.showWindow = 0;
            net = train(net, reduced_data');
            som_clusters_all(run, :) = vec2ind(net(reduced_data'));

            % Calcola il Mean Quantization Error (MQE)
            weights = net.IW{1}; % Ottieni i pesi della mappa SOM
            % Calcola le distanze tra ogni input e il rispettivo BMU
            errors = zeros(1, size(reduced_data, 1));
            for i = 1:size(reduced_data, 1)
                input_vector = reduced_data(i, :);
                distances = vecnorm(weights - input_vector, 2, 2);
                [~, bmu_index] = min(distances);
                errors(i) = distances(bmu_index);
            end
            % Calcola il MQE come la media delle distanze
            mean_quantization_error{dim_idx, epoch_idx}(run) = mean(errors);

        end

        % Calcolare AMI e Adjusted Rand Index per tutte le coppie di clusterizzazioni
        ami_values = zeros(num_runs);
        adjusted_rand_values = zeros(num_runs);

        for i = 1:num_runs
            for j = 1:num_runs
                if i < j % Salvare solo il triangolo inferiore
                    ami_values(i, j) = ami(som_clusters_all(i, :), som_clusters_all(j, :));
                    adjusted_rand_values(i, j) = rand_index(som_clusters_all(i, :), som_clusters_all(j, :), 'adjusted');
                else
                    ami_values(i, j) = NaN; % Ignorare il triangolo superiore e la diagonale
                    adjusted_rand_values(i, j) = NaN; % Ignorare il triangolo superiore e la diagonale
                end
            end
        end

        % Salvare i risultati nei cell array
        results_ami{dim_idx, epoch_idx} = ami_values(~isnan(ami_values));
        results_ari{dim_idx, epoch_idx} = adjusted_rand_values(~isnan(adjusted_rand_values));

        % Calcolare le medie dei valori di AMI e ARI
        results_mean_ami(dim_idx,epoch_idx) = nanmean(ami_values(:));
        results_mean_ari(dim_idx,epoch_idx) = nanmean(adjusted_rand_values(:));
        results_mean_qe(dim_idx,epoch_idx)  = nanmean(mean_quantization_error{dim_idx, epoch_idx});

        fprintf('\n SOM Dim: [%d %d], Epochs: %d\n', som_dim(1), som_dim(2), epochs);
        fprintf('Mean Adjusted Mutual Information (AMI): %.4f\n', results_mean_ami(dim_idx,epoch_idx));
        fprintf('Mean Adjusted Rand Index (ARI): %.4f\n', results_mean_ari(dim_idx,epoch_idx));
        fprintf('Mean Quantization Error (QE): %.4f\n', results_mean_qe(dim_idx,epoch_idx));

    end
end

%% Trova i migliori parametri per ogni som_dim
clear best_id best_epochs
for dim_idx = 1:size(som_dims, 1)
disp(' ')
    [ ~ , best_id  ] = max(results_mean_ari(dim_idx,:)) ;
    best_epochs_ari(dim_idx) = epoch_values(best_id) ;
    [ ~ , best_id  ] = max(results_mean_ami(dim_idx,:)) ;
    best_epochs_ami(dim_idx) = epoch_values(best_id) ;
    [ ~ , best_id  ] = max(results_mean_qe(dim_idx,:)) ;
    best_epochs_qe(dim_idx) = epoch_values(best_id) ;
    fprintf('For SOM Dim: [%d %d] Best Epochs (ARI) = %d\n', som_dims(dim_idx,1), som_dims(dim_idx,2), best_epochs_ari(dim_idx) );
    fprintf('For SOM Dim: [%d %d] Best Epochs (AMI) = %d\n', som_dims(dim_idx,1), som_dims(dim_idx,2), best_epochs_ami(dim_idx) );
    fprintf('For SOM Dim: [%d %d] Best Epochs (QE)  = %d\n', som_dims(dim_idx,1), som_dims(dim_idx,2), best_epochs_qe(dim_idx) );
end
disp(' ')

% Salvare i risultati per i parametri migliori
clear PipelineB_Results
PipelineB_Results.qe = mean_quantization_error ;
PipelineB_Results.ari = results_ari ;
PipelineB_Results.ami = results_ami ;
PipelineB_Results.UseDims = som_dims ;
PipelineB_Results.UseEpochs = best_epochs_ari  ;
save(['PipeB_Id' num2str(IdProcedure) '_Out_Parametrization_Ncycles' num2str(num_runs) '.mat'] , 'PipelineB_Results' )

%% Visualizzare AMI e ARI grezzi per ogni combinazione di parametri con linee e shaded error
figure('Position', [100, 100, 1200, 600]);
fsizeB = 16 ;

% Codifica dei colori per som_dims
% colors = lines(size(som_dims, 1));
colors = othercolor('Set19',size(som_dims, 1));

% Calcolare medie e errori standard per AMI e ARI
mean_ami = zeros(size(som_dims, 1), length(epoch_values));
std_ami = zeros(size(som_dims, 1), length(epoch_values));
mean_ari = zeros(size(som_dims, 1), length(epoch_values));
std_ari = zeros(size(som_dims, 1), length(epoch_values));

for dim_idx = 1:size(som_dims, 1)
    for epoch_idx = 1:length(epoch_values)
        ami = results_ami{dim_idx, epoch_idx};
        ari = results_ari{dim_idx, epoch_idx};

        mean_ami(dim_idx, epoch_idx) = mean(ami);
        std_ami(dim_idx, epoch_idx) = std(ami);
        mean_ari(dim_idx, epoch_idx) = mean(ari);
        std_ari(dim_idx, epoch_idx) = std(ari);
    end
end

% Plot AMI con linee e shaded error
subplot(1, 2, 1);
hold on;
for dim_idx = 1:size(som_dims, 1)
    epochs = epoch_values;
    ami_means = mean_ami(dim_idx, :);
    ami_stds = std_ami(dim_idx, :);

    % Plot della linea
    pl(dim_idx) = plot(epochs, ami_means, 'Color', colors(dim_idx, :), 'LineWidth', 2);

    % Plot dell'errore standard
    fill([epochs, fliplr(epochs)], [ami_means - ami_stds, fliplr(ami_means + ami_stds)], colors(dim_idx, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
xlabel('Epochs');
ylabel('Adjusted Mutual Information (AMI)');
title('Adjusted Mutual Information (AMI)');
xlim([min(epoch_values)-10 max(epoch_values)+10]);
set(gca, 'XTick', epoch_values , 'FontSize',fsizeB);
legend_entries = arrayfun(@(dim) sprintf('SOM Dim: [%d %d]', som_dims(dim, 1), som_dims(dim, 2)), 1:size(som_dims, 1), 'UniformOutput', false);
legend(pl, legend_entries, 'Location', 'southeast');

% Plot ARI con linee e shaded error
subplot(1, 2, 2);
hold on;
for dim_idx = 1:size(som_dims, 1)
    epochs = epoch_values;
    ari_means = mean_ari(dim_idx, :);
    ari_stds = std_ari(dim_idx, :);

    % Plot della linea
    pl(dim_idx) = plot(epochs, ari_means, 'Color', colors(dim_idx, :), 'LineWidth', 2);

%     Plot dell'errore standard
    fill([epochs, fliplr(epochs)], [ari_means - ari_stds, fliplr(ari_means + ari_stds)], colors(dim_idx, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
xlabel('Epochs');
ylabel('Adjusted Rand Index (ARI)');
title('Adjusted Rand Index (ARI)');
xlim([min(epoch_values)-10 max(epoch_values)+10]);
set(gca, 'XTick', epoch_values, 'FontSize',fsizeB);
legend(pl, legend_entries, 'Location', 'southeast');

sgtitle('AMI e ARI per diverse combinazioni di SOM Dim e Epochs');

flnm = ['Fig_PipeB_Id' num2str(IdProcedure) '_Out_Parametrization_Ncycles' num2str(num_runs) '_ParametersResults'] ;
export_fig(flnm, '-jpg' , '-tif') 
