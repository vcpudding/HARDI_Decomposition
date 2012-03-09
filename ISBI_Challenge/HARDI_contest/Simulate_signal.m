
clear all
close all

%% (STEP 1) Create the SAMPLING SCHEME
%  ===================================

%
% write here your own code to create the text file (i.e. "gradient_list.txt")
% containing the SAMPLING SCHEME needed by your recostruction method
%


%% (STEP 2) Simulate the SIGNAL for each voxel, phantom and SNR
%  ============================================================
SNRs          = [ 10 20 30 40 ];
PHANTOM_names = { 'Training_IV', 'Training_SF', 'Training_3D_SF' }; % in the folder "Phantoms"

% load the SAMPLING SCHEME
XYZB = dlmread( 'gradient_list.txt', ' ' );
nSAMPLES = size( XYZB, 1 );

for PHANTOM_name = PHANTOM_names
    load( fullfile('Phantoms',PHANTOM_name{1}) ) % load the variable "FIELD"
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
        save( sprintf('%s__SNR=%02d__SIGNAL.mat', PHANTOM_name{1}, SNR), 'E' )
    end
end

%
% NB: in the case of the "testing data", simply send the file "gradient_list.txt"
% to the organizers, and you'll receive back the SIGNAL simulated as in the code
% snipped above (this is the code we actually use)
%


% %% (STEP 3) Reconstruct each voxel with your own method
% %  ====================================================
% for PHANTOM_name = PHANTOM_names
% for SNR = SNRs
% 
%     % load the SIGNAL simulated with this SNR
%     load( sprintf('%s__SNR=%02d__SIGNAL.mat', PHANTOM_name{1}, SNR) )
%     [n1,n2,n3,~] = size( E );
% 
%     % reconstruction of the configuration in each voxel
%     myESTIMATION = [];
%     myESTIMATION.FIELD = cell(  n1, n2, n3 );
%     myESTIMATION.ODF   = zeros( n1, n2, n3, 724 );
% 
%     %
%     % write here your own code to apply your RECONSTRUCTION METHOD
%     % to this data and store the results in the structure "myESTIMATION"
%     % (what follows is just an illustrative example)
%     %
% 
%     for x = 1:n1
%     for y = 1:n2
%     for z = 1:n3
%         % create the MultiTensor object
%         myESTIMATION.FIELD{x,y,z} = MultiTensor();
% 
%         % number of fiber compartments in this voxel
%         myESTIMATION.FIELD{x,y,z}.M = ?????;
% 
%         % relative volume fractions and diffusivities for each estimated fiber compartment
%         myESTIMATION.FIELD{x,y,z}.f = [ ??, ??, ... ];
%         myESTIMATION.FIELD{x,y,z}.lambda = [ ?? ?? ?? ; ?? ?? ??; ... ]' * 1e-3;
% 
%         % orientation of each fiber compartment
%         myESTIMATION.FIELD{x,y,z}.R = zeros( 3, 3, myESTIMATION.FIELD{x,y,z}.M );
%         for d = 1:myESTIMATION.FIELD{x,y,z}.M
%             myESTIMATION.FIELD{x,y,z}.R(:,:,d) = myESTIMATION.FIELD{x,y,z}.ROTATION( ??, ?? );
%         end
% 
%         % estimate the ODF in this voxel
%         % (please use the directions in the file ODF_XYZ.mat)
%         myESTIMATION.ODF(x,y,z,:) = [ ????? ];
%     end
%     end
%     end
% 
%     % save the results
%     save( sprintf('%s__SNR=%02d__myESTIMATION.mat', PHANTOM_name{1}, SNR), 'myESTIMATION' );
% end
% end
% 
% %
% % NB: in the case of the "testing data", simply send these files
% % to the organizers for the evaluation of the reconstructions
% %