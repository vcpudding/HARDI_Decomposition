function testDeconvolveFibers
% close all
clear all

b=3000; % for synthetic data simulations
S0=1; 
t=pi/3; % the angle between fibers for simulations
snr=30;

load Brain_GradientOrientations
UnitVectors

% [S trueDirs]=simulateDWData(b, GradientOrientations, [0, t], [0.5,0.5], 1);
[S trueDirs]=simulateDWData(b, GradientOrientations, 0, 1, 1);

% fiber_orientation1=[cos(pi/3) sin(pi/3) 0];
% orientation1=atan2(fiber_orientation1(2),fiber_orientation1(1));
% S = Simulate_DW_data(b, GradientOrientations, orientation1, t, 1, 0);
% trueDirs = fiber_orientation1';

y=randn(length(S), 2);
S_noisy = abs(S+1/snr*(y(:,1)+sqrt(-1)*y(:,2)));
% S_noisy = S;

options.order=40; % polynomial order
options.delta= 5; %29.2; % estimated from the data, b-value multiplied by the dominant diffusivity
options.lambda=1e-6; % initial damping parameter value
options.tol=1e-6; % convergence criteria
options.maxiter=20000;
options.step = 1e-3;
options.accurate_integration = 1;

%options added by stella
options.linesearch = 0; % use line search: 1; otherwise: 0
options.innerconvergence = 0; % use inner convergence: 1; otherwise: 0
options.innertol=1e-3; % inner convergence tolerance
options.maxinneriter=1000; % max no. of inner iteration
options.init=0; % 0. random initialization; 1. DTI initialization; 2. Fixed value

disp(['snr: ', num2str(snr)]);
disp(['angle: ', num2str(t*180/pi)]);
disp(['delta: ', num2str(options.delta )]);
tic
[dirs weights stat] = deconvolveFibersGD(S_noisy, GradientOrientations', b, 3, options);
toc

disp('fiber directions');
disp(trueDirs);
disp('-------------------');
disp(dirs);
disp('-------------------');
disp(directionDeviation(dirs, trueDirs));
disp(['weights: ', num2str(weights)]);
figure(1), plot(stat.eCurve);
% trueDirs = [1 0 0; cos(t) sin(t) 0]';
% disp(['direction deviation: ', num2str(directionDeviation(dirs, trueDirs))]);

% fid = fopen(sprintf('../Data/results/log[linesearch=%d][innerconvergence=%d][snr=%d][sepAngle=%d].txt', options.linesearch, options.innerconvergence, snr, t*180/pi), 'w');
% for i=1:10
%     [dirs weights stat] = deconvolveFibers(S_noisy, GradientOrientations', b, 2, options);
%     fprintf(fid, 'Init dirs: (%0.2f, %0.2f, %0.2f; %0.2f, %0.2f, %0.2f)\n', stat.initDirs(:));
%     fprintf(fid, 'Init energy: %0.2f\n', stat.initEnergy);
%     fprintf(fid, 'Converge energy: %0.2f\n', stat.convEnergy);
%     fprintf(fid, 'Convergence : %d iterations, time = %2f\n', stat.nIts, stat.time);
%     fprintf(fid, 'Output dirs: (%0.2f, %0.2f, %0.2f; %0.2f, %0.2f, %0.2f)\n', dirs(:));
%     fprintf(fid, 'Output weights: (%0.2f, %0.2f)\n', weights);
%     fprintf(fid, '\n\n');
% end
% fclose(fid);

% deltas = [2 3 8 15 29 40];
% figure;
% for i=1:length(deltas)
%     options.delta = deltas(i);
%     [dirs weights stat] = deconvolveFibers(S_noisy, GradientOrientations', b, 2, options);
%     subplot(2, 3, i);
%     plot(stat.eCurve);
%     title(sprintf('delta = %d', deltas(i)));
%     disp(dirs);
%     disp(weights);
% end
end