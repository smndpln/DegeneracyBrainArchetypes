
fSiz = 15 ;
% ModularBeh.factors contains behavioral factors preprocessed after boxcox correction
% Here, we also normalize the factors between 1 and 100
ModularBeh.factorsNorm = normalize ( ModularBeh.factors , 'range' , [1 100] ) ;
ModularBeh.factorsLong = repmat ( ModularBeh.factorsNorm, NumRest , 1 ) ;
NumFact = size(ModularBeh.factors,2) ;
DoFig = 1 ;
DisplayAnova = 0 ;
% SetNewLabels = { 'Agreeableness' 'Openness' 'Conscentiousness' 'Neuroticism' 'Extraversion' ...
%     'DDisc'}

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

        %% ANOVA
        % Eseguire l'ANOVA a due vie con interazione
        [~, anovaTable, stats] = anovan(behav_score, {x_coords, y_coords}, 'Model', 'interaction', 'VarNames', {'X', 'Y'} , 'display',DisplayAnova);

        % Generare i test di confronti multipli : X
        [resultsX,mX,~,gnamesX] = multcompare(stats, 'Dimension', 1 ,'Display','off');
        tblX = array2table(resultsX,"VariableNames", ["Group1","Group2","Lower Limit","Difference","Upper Limit","P-value"]);
        tblX.("Group1") = gnamesX(tblX.("Group1")); tblX.("Group2") = gnamesX(tblX.("Group2"));
        tblXmeans = array2table(mX,"RowNames",gnamesX, "VariableNames",["Mean","Standard Error"]) ;
        % Generare i test di confronti multipli : Y
        [resultsY,mY,~,gnamesY] = multcompare(stats, 'Dimension', 2 ,'Display','off');
        tblY = array2table(resultsY,"VariableNames", ["Group1","Group2","Lower Limit","Difference","Upper Limit","P-value"]);
        tblY.("Group1") = gnamesY(tblY.("Group1")); tblY.("Group2") = gnamesY(tblY.("Group2"));
        tblYmeans = array2table(mY,"RowNames",gnamesY, "VariableNames",["Mean","Standard Error"]) ;
        % Generare i test di confronti multipli : XY
        [resultsXY,mXY,~,gnamesXY] = multcompare(stats, 'Dimension', [1 2] ,'Display','off');
        tblXY = array2table(resultsXY,"VariableNames", ["Group1","Group2","Lower Limit","Difference","Upper Limit","P-value"]);
        tblXY.("Group1") = gnamesXY(tblXY.("Group1")); tblXY.("Group2") = gnamesXY(tblXY.("Group2"));
        tblXYmeans = array2table(mXY,"RowNames",gnamesXY, "VariableNames",["Mean","Standard Error"]) ;

        %% per l'interazione, filtrare i confronti non interessanti
        numRows = height(tblXY);
        keepRows = true(numRows, 1);
        for i = 1:numRows
            group1 = tblXY.Group1{i};
            group2 = tblXY.Group2{i};
            % Dividere i gruppi in fattori
            factors1 = split(group1, ',');
            factors2 = split(group2, ',');
            % Contare il numero di differenze tra i fattori
            numDifferences = sum(~strcmp(factors1, factors2));
            % Se entrambi i fattori cambiano, segna la riga come non interessante
            if numDifferences == 2
                keepRows(i) = false;
            end
        end
        % Creare una nuova tabella con solo i confronti interessanti
        tblXY_interesting = tblXY(keepRows, :);
        resultsXY_interesting = resultsXY(keepRows, :);

        %%
        % Colori per SOM-Y
        colors = lines(7);
        colors = colors(3:max(x_coords)+3,:) ;

        if DoFig == 1

            % Visualizzazione dei risultati
            figure('Position', [30, 30, 450, 560]);

            % Effetto di SOM-X
            subplot(2,2,1)
            bar( tblXmeans.Mean);
            hold on;
            errorbar( tblXmeans.Mean, tblXmeans.("Standard Error"), 'k', 'LineStyle', 'none');
            ylabel([behav_label]);
            title(['X']);
            % Show Significance
            mcCell = arrayfun(@(i) resultsX(i, 1:2), 1:size(resultsX, 1), 'UniformOutput', false);
            sigstar_color(mcCell, resultsX(:,6));
            set(gca,'XTick',1:som_dim(1), 'XTickLabels',tblXmeans.Row,'YGrid','on' , 'TickLength', [0 0]);
            yl=ylim; ylim([30 yl(2)]); box off;
            set(gca,'FontSize', fSiz) ;

            % Effetto di SOM-Y
            subplot(2,2,2)
            bar( tblYmeans.Mean, 'FaceColor', 'flat');
            hold on;
            errorbar(tblYmeans.Mean, tblYmeans.("Standard Error"), 'k', 'LineStyle', 'none');
%             ylabel([behav_label]);
            title(['Y']);
            % Show Significance
            mcCell = arrayfun(@(i) resultsY(i, 1:2), 1:size(resultsY, 1), 'UniformOutput', false);
            sigstar_color(mcCell, resultsY(:,6));
            set(gca,'XTick',1:som_dim(1), 'XTickLabels',tblYmeans.Row,'YGrid','on' , 'TickLength', [0 0]);
            yl=ylim; ylim([30 yl(2)]); box off;
            set(gca,'FontSize', fSiz) ;

            % Interazione (X:Y)
            subplot(2,2,3:4)
            barGroups = reshape(tblXYmeans.Mean, som_dim);
            errorbarGroups = reshape(tblXYmeans.("Standard Error"), som_dim);
            barColors = repelem(colors, 1, 3)';
            b = bar(barGroups, 'grouped');
            for k = 1:size(barGroups, 2)
                b(k).FaceColor = colors(k, :);
            end
            hold on;
            groupCenters = nan(size(barGroups));
            for i = 1:size(barGroups, 1)
                for j = 1:size(barGroups, 2)
                    groupCenters(i, j) = b(j).XEndPoints(i);
                end
            end
            errorbar(groupCenters, barGroups, errorbarGroups, 'k', 'LineStyle', 'none');
            ylabel([behav_label] , 'Interpreter','none' , 'FontWeight','bold' );
%             title('Interaction');
            % Show Significance
            mapResults = resultsXY_interesting
            for i = 1 : som_dim(1)^2
                mapResults(resultsXY_interesting(:,1)==i,1) = groupCenters(i) ;
                mapResults(resultsXY_interesting(:,2)==i,2) = groupCenters(i) ;
            end
            mcCell = arrayfun(@(i) mapResults(i, 1:2), 1:size(mapResults, 1), 'UniformOutput', false);
            sigstar_color(mcCell, mapResults(:,6));
            set(gca,'XTick',1:3, 'XTickLabels',unique(extractBefore(tblXYmeans.Row, ",")),'YGrid','on' , 'TickLength', [0 0]);
            yl=ylim; ylim([30 yl(2)]); box off;
            legend(arrayfun(@(y) sprintf('Y=%d', y), 1:max(x_coords), 'UniformOutput', false) , 'Location','eastoutside');
            set(gca,'FontSize', fSiz) ;
            flnm = ['Fig_pipeD_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_MultComp_fact' num2str(select_fact)] ;
            export_fig(flnm, '-jpg' , '-tif', '-r500')
            
            close

        end % Eventually do Fig

        %% save
        anvs.tables{select_fact} = anovaTable ;
        anvs.stats{select_fact}  = stats ;
        anvs.Xmeans{select_fact}  = tblXmeans ;
        anvs.Ymeans{select_fact}  = tblYmeans ;
        anvs.XYmeans{select_fact} = tblXYmeans ;
        anvs.X_results{select_fact}  = tblX ;
        anvs.Y_results{select_fact}  = tblY ;
        anvs.XY_results{select_fact} = tblXY_interesting ;
        anvs.p.X(select_fact,1)  = anovaTable{2,7} ;
        anvs.p.Y(select_fact,1)  = anovaTable{3,7} ;
        anvs.p.XY(select_fact,1) = anovaTable{4,7} ;

    end

    save(['PipeD_Id' num2str(IdProcedure) '_Out_ANOVAS_Ncycles' num2str(num_runs) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '.mat'] , 'anvs' )

    %% Visualizzazione 2.0

    % % % FDR
    p_ALL = [ anvs.p.X anvs.p.Y anvs.p.XY ] ;
    % Fattore X
    pvals_X = anvs.p.X;
    [pthr_X, pcor_X, padj_X] = fdr2(pvals_X);
    padj_X_full = padj_X;
    % Fattore Y
    pvals_Y = anvs.p.Y;
    [pthr_Y, pcor_Y, padj_Y] = fdr2(pvals_Y);
    padj_Y_full = padj_Y;
    % Interazione XY
    pvals_XY = anvs.p.XY;
    [pthr_XY, pcor_XY, padj_XY] = fdr2(pvals_XY);
    padj_XY_full = padj_XY;
    % Concatena
    padj_ALL_fdr = [ padj_X_full padj_Y_full padj_XY_full ] ;

    % Bonferroni
    padj_X_bonf = anvs.p.X.*NumFact;
    padj_Y_bonf = anvs.p.Y.*NumFact;
    padj_XY_bonf = anvs.p.XY.*NumFact;
    padj_ALL_bonf = [ padj_X_bonf padj_Y_bonf padj_XY_bonf ] ;

    % Creazione delle heatmap
    figure('Position',[10 10 500 1000]);

    % Heatmap
    imagesc(p_ALL); box off;
    clim([0 1]); 
    cb = colorbar;
    cb.Label.String = 'p Value' ;
    cb.FontSize = 14 ;
    colormap("hot");
    title(['brain-behavior coding [' num2str(som_dim(1)) 'x' num2str(som_dim(1)) ']' ]);
    set(gca, 'XTick', 1:3, 'YTick', 1:NumFact , 'XTickLabel', {'X' 'Y' 'X:Y'} , 'YTickLabel', ModularBeh.factor_labs ,'TickLabelInterpreter' ,'none');
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
    flnm = ['Fig_pipeD_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_MultComp_pHeatmaps' ] ;
    export_fig(flnm, '-jpg' , '-tif', '-r500')

    %% Visualizzazione 2.0 Subselextion
    subIndices = [ 1 10 16 18 19 ...
        22 25 26 27 28 29 30 ... 
        35 36 ...
        37 38 39 40 ...
        41 42 43 ] ;
   
    % SubSelect
    p_ALL_sub = p_ALL(subIndices,:) ;
    padj_ALL_bonf_sub = padj_ALL_bonf(subIndices,:) ;
    padj_ALL_fdr_sub  = padj_ALL_fdr(subIndices,:) ;
    subLabs = {'Agreeableness' 'DelayDisc 1yr(200)' 'DelayDisc 1yr(40K)'  'DelayDisc AUC(200)' 'DelayDisc AUC(40K)' ...
        'Anger' 'Sadness' 'Life Satisfaction' 'Mean Purpose' 'Positive Affect' 'Friendship' 'Loneliness' ...
        'Perceived Stress' 'Self-Efficacy' ...
        'Cognitive: Fluid ' 'Cognitive: Early' 'Cognitive: Crystallized' 'Cognitive: Total' ...
        'Endurance' 'Gait Speed' 'Dexterity' } ;

    % Creazione delle heatmap
    figure('Position',[10 10 500 1000]);

    % Heatmap
    imagesc(p_ALL_sub); box off;
    clim([0 1]); 
    cb = colorbar;
    cb.Label.String = 'p Value' ;
    cb.FontSize = 14 ;
    colormap("hot");
    title(['brain-behavior coding [' num2str(som_dim(1)) 'x' num2str(som_dim(1)) ']' ]);
    set(gca, 'XTick', 1:3, 'YTick', 1:numel(subLabs) , 'XTickLabel', {'X' 'Y' 'X:Y'} , 'YTickLabel', subLabs ,'TickLabelInterpreter' ,'none');
    hold on;

    % Aggiungi simboli "*" e linee bianche
    xl=xlim; yl=ylim;
    for i = 1:numel(subLabs)
        for j = 1:3
            if padj_ALL_bonf_sub(i, j) < 0.05
                plot(j,i,'Marker','o' , 'MarkerFaceColor','y','MarkerEdgeColor','y','MarkerSize',6); %(j, i, '*', 'HorizontalAlignment', 'Center', 'Color', 'y' , 'FontSize',fSiz+7);
            elseif padj_ALL_fdr_sub(i, j) < 0.05
                plot(j,i,'Marker','o' , 'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5); %(j, i, '*', 'HorizontalAlignment', 'Center', 'Color', 'y' , 'FontSize',fSiz+7);
            elseif p_ALL_sub(i, j) < 0.05
                plot(j,i,'Marker','o' , 'MarkerFaceColor',[.6 .6 .6],'MarkerEdgeColor',[.6 .6 .6],'MarkerSize',3); %(j, i, '*', 'HorizontalAlignment', 'Center', 'Color', 'y' , 'FontSize',fSiz+7);
            end
            plot(xl,[i i]-.5,'Color','w' , 'LineWidth',1)
            plot([j j]+.5,yl,'Color','w' , 'LineWidth',1)
        end
    end
    box off;
    set(gca,'FontSize',fSiz-1)

    % Salva le heatmap come immagine
    flnm = ['Fig_pipeDsub_Id' num2str(IdProcedure) '_Dim' num2str(som_dim(1)) '_epochs' num2str(som_epochs) '_MultComp_pHeatmaps' ] ;
    export_fig(flnm, '-jpg' , '-tif', '-r500')

end

%% RUN analysis on latent factors ?
% Pipeline_D_2wayANOVA_latent

