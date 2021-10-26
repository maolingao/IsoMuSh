% demo: tosca centaur
%
clear; restoredefaultpath
addpath(genpath('./code'))
addpath(genpath('./eval'))


%% choose shape collection
% For TOSCA centaur we provide the original *.off and the corresponding
% preprocessed *.mat files used in our pipeline for demo purpose.
% To experiment with other shape collections, please download the original
% TOSCA dataset, and adapt params.dataPreparation.datapath_orig.

datasetname = 'tosca';
dataname = 'centaur';
params.dataPreparation.datapath_orig = './data/tosca/off_orig'; %'<ROOT_TO_TOSCA_OFF_FILES>';


%% safe to ignore

params.verbose = 2;
params.isomush.maxIter = 1e2;
numFace = 2e3;
dimLB = 1e2;

%% do not touch
params.rootpath = pwd;
params.datasetname = datasetname;
params.dataname = dataname;
params.dataPreparation.dimLB = dimLB;
params.dataPreparation.numFace = numFace;
params.dataPreparation.datapath_mat = fullfile(params.rootpath, 'data', datasetname, 'mat_orig');

%%
% load dataset
fprintf('loading (%s, %s)... \n', datasetname, dataname);
data = load_dataset(params);

% initialisation
% -- ZoomOut --
fprintf('\nInitialisation...\n')
fprintf('\n1) Running ZoomOut for all pairs on (%s, %s)... \n', datasetname, dataname);
[C, C_map, time_zoomout] = zoomout_wrapper(data, params);


% -- Synchronisation --
fprintf('\n2) Synchronising pairwise results on (%s, %s)... \n', datasetname, dataname);
[U, Q] = synchronisation(C, data, params);



% IsoMush
fprintf('\nStarting IsoMuSh...\n')
fprintf('\nRunning IsoMuSh on (%s, %s)... \n', datasetname, dataname);
[U, Q] = isomush(U, Q, data, params);








