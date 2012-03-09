
close all;
clear all;

addpath('../Program');

fiber_orientation1=[0 1 0];
orientation1=atan2(fiber_orientation1(2),fiber_orientation1(1));
 
b=2000; % for synthetic data simulations
S0=1; 
t=pi/2; % the angle between fibers for simulations
SNR=100;
sig=1/SNR; % SNR for synetheric data

options.order=16; % polynomial order
options.delta=3; % estimated from the data, b-value multiplied by the dominant diffusivity
options.lambda=1e-3; % initial damping parameter value
options.tol=1e-6; % convergence criteria
options.maxiter=20000;

%%% these parameters are not in use %%%%%%%%%%%%%%%%%%%%%%%

options.rank=2; % initial rank (number of fibers)
options.modsel=1e-10; % thereshold for model selection
options.init=0; % random initialization, 1 for signal minima initialization
multiple_init=1; % multiple random initialization
init_number=10; % number of multiple starts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

accurate_integration=0;

if accurate_integration
    load sphere_integration_2562
else 
    load sphere_integration_642
end

g2=centroids; 

load Brain_GradientOrientations

% The convolution kernel - delta is estimated from the data
kernel=exp(-options.delta*(GradientOrientations*g2').^2);

  
% synthetic data simulation 
R=[cos(-t) sin(-t) 0;-sin(-t) cos(-t) 0;0 0 1];
fiber_orientation2=fiber_orientation1*R';
% [S trueDirs]=simulateDWData(b, GradientOrientations, [0 t], [0.5, 0.5], [1.7e-3, 0.3e-3, 0.3e-3]);
% [S trueDirs]=simulateDWData(b, GradientOrientations, 0, 1, [1.7e-3, 0.3e-3, 0.3e-3]);
S = Simulate_DW_data(b, GradientOrientations, orientation1, t, 1, 0);

y=randn(2,length(S));
S_noisy = abs(S'+sig*(y(1,:)+sqrt(-1)*y(2,:)));

Ev=randn(3,3);
Ev=Stabilized_Gram_Schmidt(Ev); % orthonormalize

disp(['snr: ', num2str(SNR)]);
disp(['angle: ', num2str(t*180/pi)]);
disp(['delta: ', num2str(options.delta )]);

[weights fiber_directions final_func]=deconvolve_two_fibers(S_noisy, kernel, centroids, tri_areas, GradientOrientations, Ev(:,1:2), options);
disp(fiber_directions);
disp(weights);
