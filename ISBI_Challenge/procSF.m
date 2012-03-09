
close all;
clear all;
addpath ../../HARDI_Denoising/Program/HARDI_simulations
addpath ../Program
addpath HARDI_contest
load Brain_GradientOrientations

SNRs          = [ 10 20 30 ];
PHANTOM_names = { 'Training_SF', 'Training_3D_SF' }; % in the folder "Phantoms"
bVal          = 2000;
order         = 16;

for PHANTOM_name = PHANTOM_names
for SNR = SNRs
for denoise = 0

    load( fullfile('TrainingData/Phantoms',PHANTOM_name{1}) );
    trueField = FIELD;
    
    % load the SIGNAL simulated with this SNR
    load( sprintf('TrainingData/%s__B=%04d__SNR=%02d__SIGNAL.mat', PHANTOM_name{1}, bVal, SNR) );
    [n1,n2,n3,nGrads] = size( E );
    % anisotropic denoising
    if denoise
        E = dwiAnisotropicFiltering(double(E), 1/SNR, 0.5, 3000, 1e-3, 1e-6);
    else
        E = double(E);
    end

    % reconstruction of the configuration in each voxel
    myESTIMATION = [];
    myESTIMATION.FIELD = cell(  n1, n2, n3 );
    myESTIMATION.ODF   = zeros( n1, n2, n3, 724 );

    %
    % write here your own code to apply your RECONSTRUCTION METHOD
    % to this data and store the results in the structure "myESTIMATION"
    % (what follows is just an illustrative example)
    %
    
    logFileID = fopen(sprintf('Results/%s__B=%04d__SNR=%02d__DELTA=auto__DENOISE=%d_GD.txt', PHANTOM_name{1}, bVal, SNR, denoise), 'w');

    for x = 1:n1
    for y = 1:n2
    for z = 1:n3
        
        S = double(E(x,y,z,:));
        S = reshape(S, [numel(S), 1]);
        deconvOptions.order=order; % polynomial order
        deconvOptions.lambda=1e-3; % initial damping parameter value
        deconvOptions.tol=1e-5; % convergence criteria
        deconvOptions.maxiter=8000;
        deconvOptions.step = 1e-3;
        deconvOptions.accurate_integration = 0;

        %options added by stella
        deconvOptions.linesearch = 1; % use line search: 1; otherwise: 0
        deconvOptions.innerconvergence = 0; % use inner convergence: 1; otherwise: 0
        deconvOptions.innertol=1e-3; % inner convergence tolerance
        deconvOptions.maxinneriter=100; % max no. of inner iteration
        deconvOptions.init=0; % 0. random initialization; 1. DTI initialization
        
        deconvOptions.delta = 10;
        deltaStep = 10;

%         [dirs, weights, stat] = deconvolveFibersGD(S, GradientOrientations', bVal, 3, deconvOptions);
        
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
        
        for it=1:10
            disp(['delta: ', num2str(deconvOptions.delta)]);
            [dirs, weights, stat] = estimateFibers(S, GradientOrientations', bVal, estimateOptions, deconvOptions);
            if ~isempty(weights)
                break;
            end
            deconvOptions.delta = deconvOptions.delta+deltaStep;
        end
        
        if isempty(weights)
            continue;
        end
        
        % create the MultiTensor object
        myESTIMATION.FIELD{x,y,z} = MultiTensor();

        % number of fiber compartments in this voxel
        myESTIMATION.FIELD{x,y,z}.M = length(weights);

        % relative volume fractions and diffusivities for each estimated fiber compartment
        myESTIMATION.FIELD{x,y,z}.f = weights;
        myESTIMATION.FIELD{x,y,z}.lambda = [ 0.3 0.3 1.7 ]' * 1e-3;

        % orientation of each fiber compartment
        myESTIMATION.FIELD{x,y,z}.R = zeros( 3, 3, myESTIMATION.FIELD{x,y,z}.M );
        for d = 1:myESTIMATION.FIELD{x,y,z}.M
            myESTIMATION.FIELD{x,y,z}.R(:,:,d) = myESTIMATION.FIELD{x,y,z}.ROTATION( atan2(dirs(2,d), dirs(1,d)), acos(dirs(3,d)) );
        end        
        
        nTrueDirs = trueField{x,y,z}.M;
        trueDirs = reshape(trueField{x,y,z}.R(:,3,:), [3, nTrueDirs]);
        %dev = directionDeviation(dirs, trueDirs);
        dev = calc_deviation(dirs, trueField{x,y,z}.R);
        disp('ground truth');
        disp(trueField{x,y,z}.f);
        disp('estimated');
        disp(myESTIMATION.FIELD{x,y,z}.f);
        disp('deviation');
        disp(dev);

        % estimate the ODF in this voxel
        % (please use the directions in the file ODF_XYZ.mat)
%         myESTIMATION.ODF(x,y,z,:) = [ ????? ];
        fprintf(logFileID, '%d, ', [x, y, z]);
        fprintf(logFileID, '\t\t');
        fprintf(logFileID, '%d, \t\t%#0.4g,\t\t%#0.4g,\t\t%d-%d', deconvOptions.delta, stat.convEnergy, mean(dev), length(weights), nTrueDirs);
        fprintf(logFileID, '\n');
    end
    end
    end
    
    fclose(logFileID);

    % save the results
    save( sprintf('Results/%s__B=%04d__SNR=%02d__DELTA=auto__DENOISE=%d_GD.mat', PHANTOM_name{1}, bVal, SNR, denoise), 'myESTIMATION' );
end
end
end
