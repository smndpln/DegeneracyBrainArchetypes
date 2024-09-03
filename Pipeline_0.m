
disp('Pipeline_A_ImportData')
Pipeline_A_ImportData

% Insert ID of the procedure
IdProcedure = 99 ;

% Good Parameters for SOM Parametrization
num_runs = 100 ;
som_dims = [ 2 2 ];
epoch_values = 100 ; 
% % % % % USED FOR PRETRAINING (FINE TUNING)
% % % % % som_dims = [2 2; 3 3; 4 4; 5 5; 6 6 ; 7 7; 8 8; 9 9; 10 10];
% % % % % epoch_values = [25:25:400]; 

disp('Pipeline_B_SomParametization')
Pipeline_B_SomParametization

PipelineB_Results.UseEpochs = [ 200 200 ] ;

disp('Pipeline_C_ConsensusSOM')
Pipeline_C_ConsensusSOM

disp('Pipeline_D_2wayANOVA')
Pipeline_D_2wayANOVA

disp('Pipeline_D_MixedModel')
Pipeline_D_MixedModel



