
% Please note that these script consider the parcellation from Joliot and colleagues (2015) and the
% network stucture from Di Plinio & Ebisch, 2018

% Joliot, M., Jobard, G., Naveau, M., Delcroix, N., Petit, L., Zago, L., et al. (2015). AICHA: an atlas of intrinsic connectivity of homotopic areas. J. Neurosci. Methods 254, 46–59. doi: 10.1016/j.jneumeth.2015.07.013
% Di Plinio, S., and Ebisch, S. J. H. (2018). Brain network profiling defines functionally specialized cortical networks. Hum. Brain Mapp. 39, 4689–4706. doi: 10.1002/hbm.24315

%% A: DATA IMPORT

% PERSONALITY dataset
% load('***.mat')  % Import Brain data
% load('***.mat')  % Import Behavioral data

load('XYZ_AI.mat')
set(0,'defaultfigurecolor',[1 1 1])
rng(10) ;

% PARAMETERS
numSubjects = size(FC.z,3) ;        % Choose number of subjects
UseGrid = 'gridtop' ;       %    hextop   gridtop

UseLatent = 0 ;
NumRest = size(FC.z,4) ;
Nregions = 346 ; % using Joliot's parcellation

addpath(genpath('C:\Users\Hp\Documents\MATLAB\BCT\'))
addpath(genpath('C:\Users\Hp\Documents\MATLAB\Export_Fig\'))
addpath(genpath('C:\Users\Hp\Documents\MATLAB\othercolor\'))
addpath(genpath('C:\Users\Hp\Documents\MATLAB\subplot_tight\'))

%% BRAIN VARIABLE
connectivity_matrices = [] ;
% Prendo solo la corteccia 
for r = 1 : NumRest
    preconnectivity_matrices = FC.z(1:Nregions,1:Nregions,1:numSubjects,r);
    connectivity_matrices = cat ( 3 ,connectivity_matrices , preconnectivity_matrices) ;
end
connectivity_matrices (abs(connectivity_matrices)>4) = NaN ;

%% BNP Stuff
% Qui importiamo la suddivisione in network di Di Plinio & Ebisch (2018)
% per una più chiara interpretazione dei cluster/neuroni
load('ConversionMat.mat'); % Uso ConversionMat

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

%% Riduzione della Dimensionalità con PCA
% Estrazione e appiattimento della parte inferiore triangolare di ciascuna matrice
numLowerTriElems = Nregions*(Nregions-1)/2;
flattened_matrices = zeros(numSubjects, numLowerTriElems);

for n = 1:numSubjects*NumRest
    lowerTri = tril(connectivity_matrices(:,:,n), -1);
    flattened_matrices(n, :) = lowerTri(lowerTri ~= 0);
end

% PCA
[coeff, score, ~, ~, explained] = pca(flattened_matrices);
numComponents = 14 ; % based on parallel analysis!!! pa_test(flattened_matrices(1:100,:)) ;
reduced_data = zscore(score(:, 1:numComponents)); % QUESTI SARANNO I DATI USATI IN SOM!

