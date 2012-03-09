function result = procVoxel (phantomName, x, y, z, bVal, SNR, denoise)
        
addpath ../Program
addpath HARDI_contest
load Brain_GradientOrientations

folderName = sprintf('Results/%s__B=%04d__SNR=%02d__DELTA=auto__DENOISE=%d', phantomName, bVal, SNR, denoise);
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

logFileID = fopen(fullfile(folderName, 'parlog.txt'), 'a');
load( fullfile('TrainingData/Phantoms',phantomName) );
trueField = FIELD;

if denoise
    fileName = sprintf('TrainingData/%s__B=%04d__SNR=%02d__SIGNALfil.mat', phantomName, bVal, SNR);
    if ~exist(fileName, 'file')
        load (sprintf('TrainingData/%s__B=%04d__SNR=%02d__SIGNAL.mat', phantomName, bVal, SNR));
        E = dwiAnisotropicFiltering(double(E), 1/SNR, 0.5, 3000, 1e-6, 1e-6);
        E = single(E);
        save(fileName, 'E');
    end
    load( fileName );
else
    load( sprintf('TrainingData/%s__B=%04d__SNR=%02d__SIGNAL.mat', phantomName, bVal, SNR) );
end

S = double(E(x,y,z,:));
S = reshape(S, [numel(S), 1]);
deconvOptions.order=16; % polynomial order
deconvOptions.lambda=1e-3; % initial damping parameter value
deconvOptions.tol=1e-6; % convergence criteria
deconvOptions.maxiter=20000;
deconvOptions.step = 1e-3;
deconvOptions.accurate_integration = 0;

%options added by stella
deconvOptions.linesearch = 1; % use line search: 1; otherwise: 0
deconvOptions.innerconvergence = 0; % use inner convergence: 1; otherwise: 0
deconvOptions.innertol=1e-3; % inner convergence tolerance
deconvOptions.maxinneriter=100; % max no. of inner iteration
deconvOptions.init=0; % 0. random initialization; 1. DTI initialization
deconvOptions.delta = 10;
deconvOptions.deltaStep = 10;

switch SNR
    case 10
        estimateOptions.energy_thres = 0.95;
    case 20
        estimateOptions.energy_thres = 0.6;
    case 30
        estimateOptions.energy_thres = 0.6;
end

estimateOptions.min_num_valid = 3;
estimateOptions.min_num_trial = 3;
estimateOptions.max_num_trial = 5;
estimateOptions.algorithm = 0;
estimateOptions.fraction_thres = 1e-2;


for it=1:5
    disp(['delta: ', num2str(deconvOptions.delta)]);
    [dirs, weights, stat] = estimateFibers(S, GradientOrientations', bVal, estimateOptions, deconvOptions);
    if ~isempty(weights)
        break;
    end
    deconvOptions.delta = deconvOptions.delta+deconvOptions.deltaStep;
end

result = MultiTensor();
        
if isempty(weights)
    deconvOptions.delta = 10;
    [dirs, weights, stat] = deconvolveFibersGD(S, GradientOrientations', bVal, 3, deconvOptions);
    idx = weights>estimateOptions.fraction_thres;
    weights = weights(idx);
    weights = weights/sum(weights);
    dirs = dirs(:, idx);
end

% create the MultiTensor object

% number of fiber compartments in this voxel
result.M = length(weights);

% relative volume fractions and diffusivities for each estimated fiber compartment
result.f = weights;
result.lambda = [ 0.3 0.3 1.7 ]' * 1e-3;
result.lambda = result.lambda(:, ones(1,result.M));

% orientation of each fiber compartment
result.R = zeros( 3, 3, result.M );
for d = 1:result.M
    result.R(:,:,d) = result.ROTATION( atan2(dirs(2,d), dirs(1,d)), acos(dirs(3,d)) );
end

fileName = fullfile(folderName, sprintf('v_%02d_%02d_%02d.mat', x,y,z));
save(fileName, 'result');


nTrueDirs = trueField{x,y,z}.M;
dev = calc_deviation(dirs, trueField{x,y,z}.R);
fprintf(logFileID, '%d, ', [x, y, z]);
fprintf(logFileID, '\t\t');
fprintf(logFileID, '%d, \t\t%#0.4g,\t\t%#0.4g,\t\t%d-%d', deconvOptions.delta, stat.convEnergy, mean(dev), length(weights), nTrueDirs);
fprintf(logFileID, '\n');

end
