% Pipeline_Z_Control_AssignationsHistogram

NumClusters = som_dim(1)*som_dim(2) ;

% Simuliamo i dati dei cluster per ogni soggetto e misura ripetuta.
% I cluster sono rappresentati da numeri interi da 1 a NumClusters.
data = randi(NumClusters, numSubjects, NumRest);

% Inizializziamo un array per memorizzare il numero di cluster unici per ogni soggetto
uniqueClustersCount = zeros(numSubjects, 1);

% Per ogni soggetto, contiamo il numero di cluster unici a cui appartiene
for i = 1:numSubjects
    uniqueClustersCount(i) = length(unique(data(i, :)));
end

% Se vuoi visualizzare una tabella con i risultati
simDataTable = table((1:numSubjects)', uniqueClustersCount, ...
    'VariableNames', {'Soggetto', 'NumeroClusterUnici'});

% disp(simDataTable);
