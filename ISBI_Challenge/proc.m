addpath ../Program
load Brain_GradientOrientations

%% (STEP 1) Create the SAMPLING SCHEME
%  ===================================

%
% write here your own code to create the text file (i.e. "gradient_list.txt")
% containing the SAMPLING SCHEME needed by your recostruction method
%

XYZB = [GradientOrientations, bVal*ones(length(GradientOrientations),1)];
dlmwrite('gradient_list.txt', XYZB, ' ');


%% (STEP 2) Simulate the SIGNAL for each voxel, phantom and SNR
%  ============================================================
SNRs          = [ 10 20 30 40 ];
PHANTOM_names = { 'Training_IV', 'Training_SF', 'Training_3D_SF' }; % in the folder "Phantoms"
bVal = 2000;

% load the SAMPLING SCHEME
XYZB = dlmread( 'gradient_list.txt', ' ' );
nSAMPLES = size( XYZB, 1 );

for PHANTOM_name = PHANTOM_names
    load( fullfile('TrainingData/Phantoms',PHANTOM_name{1}) ) % load the variable "FIELD"
    [n1,n2,n3] = size( FIELD );

    for SNR = SNRs
        sigma = 1 / SNR;
        E = zeros( [n1 n2 n3 nSAMPLES], 'single' );
        for x = 1:n1
        for y = 1:n2
        for z = 1:n3
            VOXEL = FIELD{x,y,z};
            E(x,y,z,:) = VOXEL.acquireWithScheme( XYZB, sigma );
        end
        end
        end
        save( sprintf('TrainingData/%s__B=%04d__SNR=%02d__SIGNAL.mat', PHANTOM_name{1}, bVal, SNR), 'E' )
    end
end

%
% NB: in the case of the "testing data", simply send the file "gradient_list.txt"
% to the organizers, and you'll receive back the SIGNAL simulated as in the code
% snipped above (this is the code we actually use)
%


%% (STEP 3) Reconstruct each voxel with your own method
%  ====================================================
SNRs          = [ 10 20 30 40 ];
PHANTOM_names = { 'Training_IV' }; %, 'Training_SF', 'Training_3D_SF' }; % in the folder "Phantoms"
DELTAs        = [ 10 20 30 ];
nRep          = 5;
bVals         = [ 1500 3000 ];

for PHANTOM_name = PHANTOM_names
for bVal = bVals
for SNR = SNRs
for delta = DELTAs

    load( fullfile('TrainingData/Phantoms',PHANTOM_name{1}) );
    trueField = FIELD;
    
    % load the SIGNAL simulated with this SNR
    load( sprintf('TrainingData/%s__B=%04d__SNR=%02d__SIGNAL.mat', PHANTOM_name{1}, bVal, SNR) );
    [n1,n2,n3,nGrads] = size( E );
    n1 = min(n1, 10);
    n2 = min(n2, 10);
    n3 = min(n3, 10);

    % reconstruction of the configuration in each voxel
    myESTIMATION = [];
    myESTIMATION.FIELD = cell(  n1, n2, n3 );
    myESTIMATION.ODF   = zeros( n1, n2, n3, 724 );

    %
    % write here your own code to apply your RECONSTRUCTION METHOD
    % to this data and store the results in the structure "myESTIMATION"
    % (what follows is just an illustrative example)
    %
    
    logFileID = fopen(sprintf('Results/%s__B=%04d__SNR=%02d__DELTA=%02d__ESTIMATE__log.txt', PHANTOM_name{1}, bVal, SNR, delta), 'w');

    for x = 1:n1
    for y = 1:n2
    for z = 1:n3
    for iRep = 1:nRep+1
        
        S = E(x,y,z,:);
        S = reshape(S, [numel(S), 1]);
        options.order=16; % polynomial order
        options.lambda=1e-3; % initial damping parameter value
        options.tol=1e-5; % convergence criteria
        options.maxiter=8000;
        options.step = 1e-3;
        options.accurate_integration = 0;

        %options added by stella
        options.linesearch = 1; % use line search: 1; otherwise: 0
        options.innerconvergence = 1; % use inner convergence: 1; otherwise: 0
        options.innertol=1e-3; % inner convergence tolerance
        options.maxinneriter=100; % max no. of inner iteration
        if iRep<=nRep
            options.init=0; % 0. random initialization; 1. DTI initialization
        else
            options.init=1; % 0. random initialization; 1. DTI initialization
        end
        
        options.delta = delta;

        [dirs, weights, stat] = deconvolveFibersGD(double(S), GradientOrientations', bVal, 3, options);
       
%         dirs = dirs(:, weights>0.1);
%         weights = weights(weights>0.1);
        
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
        dev = directionDeviation(dirs, trueDirs);
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
        fprintf(logFileID, '%0.4f,\t\t', stat.convEnergy);
        fprintf(logFileID, '%#0.4g, %#0.4g,\t\t', mean(dev), dot(dev, weights));
        fprintf(logFileID, '% 0.4f, ', dirs(:));
        fprintf(logFileID, '\t\t');
        fprintf(logFileID, '%0.4f, ', weights(:));
        fprintf(logFileID, '\n');
    end
    end
    end
    end
    
    fclose(logFileID);

    % save the results
    save( sprintf('Results/%s__B=%04d__SNR=%02d__DELTA=%02d__myESTIMATION.mat', PHANTOM_name{1}, bVal, SNR, delta), 'myESTIMATION' );
end
end
end
end

%
% NB: in the case of the "testing data", simply send these files
% to the organizers for the evaluation of the reconstructions
%