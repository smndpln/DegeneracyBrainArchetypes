
fSiz = 14 ;
% ModularBeh.factors contains factors preprocessed after boxcox correction
% Here, we also normalize the factors between 1 and 100
ModularBeh.factorsNorm = normalize ( ModularBeh.factors , 'range' , [1 100] ) ;
ModularBeh.factorsLong = repmat ( ModularBeh.factorsNorm, NumRest , 1 ) ;
NumFact = size(ModularBeh.factors,2) ;
DoFig = 1 ;
DisplayAnova = 1 ;

%% Behav
for choose_index = 1 : 2 %size(PipelineB_Results.UseDims,1)

    som_clusters = best_clustering{choose_index} ;
    som_dim    = PipelineB_Results.UseDims(choose_index,:) ;
    som_epochs = PipelineB_Results.UseEpochs(choose_index) ;
    clear anv

    for select_fact = 1 : NumFact

        behav_score = ModularBeh.factorsLong(:,select_fact) ;
        behav_label = ModularBeh.factor_labs{select_fact} ;

        %% Generate grid
        % Define SOM dimensions and cluster coordinates for a 2x2 or 3x3 grid
        if choose_index==1
            cluster_coordinates = [
                1, 1; % Cluster 1
                2, 1; % Cluster 2
                1, 2; % Cluster 3
                2, 2; % Cluster 4
                ];
        elseif choose_index==2
            cluster_coordinates = [
                1, 1; % Cluster 1
                2, 1; % Cluster 2
                3, 1; % Cluster 3
                1, 2; % Cluster 4
                2, 2; % Cluster 5
                3, 2; % Cluster 6
                1, 3; % Cluster 7
                2, 3; % Cluster 8
                3, 3  % Cluster 9
                ];
        end

        % Map the clusters to their coordinates
        x_coords = cluster_coordinates(som_clusters, 1);
        y_coords = cluster_coordinates(som_clusters, 2);

        % Create a table with the data
        T = table(categorical(som_clusters), categorical(x_coords), categorical(y_coords), behav_score, ...
            'VariableNames', {'Cluster', 'X', 'Y', 'Behavior'});

        % Fit a linear mixed-effects model including the topological random effects
        lme1 = fitlm(T, 'Behavior~ 1 + X*Y ');
        lme1_anv = lme1.anova ;

        % Eseguire l'ANOVA a due vie con interazione
        [~, anovaTable, stats] = anovan(behav_score, {x_coords, y_coords}, 'Model', 'interaction', 'VarNames', {'X', 'Y'} ,'display','off' );

        %% save
        lmems.lme{select_fact,1}  = lme1 ;
        lmems.p.X(select_fact,1)  = lme1_anv.pValue(1) ;
        lmems.p.Y(select_fact,1)  = lme1_anv.pValue(2) ;
        lmems.p.XY(select_fact,1) = lme1_anv.pValue(3);
        lmems.anvp.X(select_fact,1)  = anovaTable{2,7} ;
        lmems.anvp.Y(select_fact,1)  = anovaTable{3,7} ;
        lmems.anvp.XY(select_fact,1) = anovaTable{4,7} ;

    end

    save(['PipeDlme_Id' num2str(IdProcedure) '_Out_ANOVAS_Ncycles' num2str(num_runs) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '.mat'] , 'lmems' )

    %% Visualizzazione 2.0

    % % % FDR
    p_ALL = [ lmems.p.X lmems.p.Y lmems.p.XY ] ;
    % Fattore X
    pvals_X = lmems.p.X;
    [pthr_X, pcor_X, padj_X] = fdr2(pvals_X);
    padj_X_full = padj_X;
    % Fattore Y
    pvals_Y = lmems.p.Y;
    [pthr_Y, pcor_Y, padj_Y] = fdr2(pvals_Y);
    padj_Y_full = padj_Y;
    % Interazione XY
    pvals_XY = lmems.p.XY;
    [pthr_XY, pcor_XY, padj_XY] = fdr2(pvals_XY);
    padj_XY_full = padj_XY;
    % Concatena
    padj_ALL_fdr = [ padj_X_full padj_Y_full padj_XY_full ] ;

    % Bonferroni
    padj_X_bonf = lmems.p.X.*NumFact;
    padj_Y_bonf = lmems.p.Y.*NumFact;
    padj_XY_bonf = lmems.p.XY.*NumFact;
    padj_ALL_bonf = [ padj_X_bonf padj_Y_bonf padj_XY_bonf ] ;

    % Creazione delle heatmap
    figure('Position',[20 20 600 800]);

    % Heatmap
    imagesc(p_ALL); box off;
    clim([0 1]);
    cb = colorbar;
    cb.Label.String = 'p Value' ;
    cb.FontSize = 14 ;
    colormap("hot");
    title(['brain-behavior coding [' num2str(som_dim(1)) 'x' num2str(som_dim(1)) ']' ]);
    set(gca, 'XTick', 1:3, 'YTick', 1:NumFact , 'XTickLabel', {'X' 'Y' 'X:Y'} , 'YTickLabel', ModularBeh.factor_labs );
    hold on;

    % Aggiungi simboli "*" e linee bianche
    xl=xlim; yl=ylim;
    for i = 1:NumFact
        for j = 1:3
            if padj_ALL_bonf(i, j) < 0.05
                plot(j,i,'Marker','o' , 'MarkerFaceColor','y','MarkerEdgeColor','y','MarkerSize',6); %(j, i, '*', 'HorizontalAlignment', 'Center', 'Color', 'y' , 'FontSize',fSiz+7);
            elseif padj_ALL_fdr(i, j) < 0.05
                plot(j,i,'Marker','o' , 'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5); %(j, i, '*', 'HorizontalAlignment', 'Center', 'Color', 'y' , 'FontSize',fSiz+7);
            elseif p_ALL(i, j) < 0.05
                plot(j,i,'Marker','o' , 'MarkerFaceColor',[.6 .6 .6],'MarkerEdgeColor',[.6 .6 .6],'MarkerSize',3); %(j, i, '*', 'HorizontalAlignment', 'Center', 'Color', 'y' , 'FontSize',fSiz+7);
            end
            plot(xl,[i i]-.5,'Color','w' , 'LineWidth',1)
            plot([j j]+.5,yl,'Color','w' , 'LineWidth',1)
        end
    end
    box off;
    set(gca,'FontSize',fSiz-1)

    % Salva le heatmap come immagine
    flnm = ['Fig_pipeDlme_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_MultComp_pHeatmaps' ] ;
    export_fig(flnm, '-jpg' , '-tif')

end

