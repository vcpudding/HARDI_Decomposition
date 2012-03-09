
clear all
close all

addpath C:\Users\Yaniv\Documents\Postdoc\MATLAB\4thorder\Decomposition_project
addpath C:\Users\Yaniv\Documents\Postdoc\MATLAB\Qball
addpath C:\Users\Yaniv\Documents\Postdoc\MATLAB\4thorder\tensor_toolbox_2.4 
addpath C:\Users\Yaniv\Documents\Postdoc\MATLAB\4thorder\tensor_toolbox_2.4\algorithms 
addpath C:\Users\Yaniv\Documents\Postdoc\MATLAB\L1magic\Optimization
addpath C:\Users\Yaniv\Documents\Postdoc\MATLAB\4thorder\Decomposition_project\decomposition_low_rank
addpath C:\Users\Yaniv\Documents\Postdoc\MATLAB\4thorder\Decomposition_project\decomposition_low_rank\main_code
addpath C:\Users\Yaniv\Documents\Postdoc\MATLAB\DTI_Estimation
addpath C:\Users\Yaniv\Documents\Postdoc\MATLAB\CSD


%addpath /scratch/Documents/Postdoc/MATLAB/CSD

global S kernel g2 tri_areas options

%%%%%%%%%%%%%%% options for low-rank spherical deconvolution %%%%%%%%%%%%
load sphere_integration_642
options.order=12; % polynomial order
options.delta=6; % estimated from the data, b-value multiplied by the dominant diffusivity
options.lambda=1e-3; % initial damping parameter value
options.tol=1e-6; % convergence criteria
options.rank=2; % initial rank (number of fibers)
options.modsel=1e-10; % thereshold for model selection
options.maxiter=200;
options.init=0; % random initialization, 1 for signal minima initialization
NewDirections=BuildSphere(1);
v=BuildSphere(6);


%S1=DWI_to_SH1(S,GradientOrientations,NewDirections,12,0.006,false,0.01)'; % convert to SH and regularize
%GradientOrientations=NewDirections;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = RandStream.create('mt19937ar','seed',5489);
RandStream.setDefaultStream(s);
defaultStream = RandStream.getDefaultStream;
savedState = defaultStream.State;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% structured field data
                                           
E0=1;
SNR=40;
order=4;
delta=200;
sigma=1/SNR;

VOXEL = MultiTensor();
myODF = MultiTensor();

load( fullfile('TrainingData','ODF_XYZ') )
%load( fullfile('TrainingData','Phantoms','Training_SF') )
%load Training_3D_SF
load Training_SF

load gradient_list_64_3000

Gradients(:,4)=1500; % change b-value to 1500

[n1 n2 n3]=size(FIELD);

g2=centroids;

%Gradients=[NewDirections, ones(length(NewDirections),1)*3000];

%%% constrined spherical deconvolution %%%%%%%%%%%%%
pscheme = gen_scheme (BuildSphere(6), 16); % was v
%DW_scheme = gen_scheme ('dir60.txt', 8);
DW_scheme = gen_scheme (Gradients(:,1:3), 8);
DW_scheme1 = gen_scheme (Gradients(:,1:3), 12);
HR_scheme = gen_scheme ('dir300.txt', 12);
HR_scheme1 = DW_scheme1;
R_SH = amp2SH (eval_DT (0.8, Gradients(1,4)/1000, pscheme), pscheme);
R_RH = SH2RH (R_SH);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kernel=exp(-options.delta*(Gradients(:,1:3)*centroids').^2); % kernel for low-rank spherical deconvolution
kernel1=exp(-options.delta*(HR_scheme.vert*centroids').^2); % kernel for low-rank spherical deconvolution

idx=1; gamma=0.01;
for k=0:2:12, % the maximal order of the HR_scheme
    v=ones(1,2*k+1);
    Lv(idx:idx+size(v,2)-1)=k*v;
    idx=idx+size(v,2);
end

L=diag(Lv.^2.*(Lv+1).^2);


for x=1:n1
    for y=1:n2 %1:n2 %n2 %2:2 %1:n2
    VOXEL=FIELD{x,y,1};
    defaultStream.State = savedState;
    E(:,x,y) = VOXEL.acquireWithScheme( Gradients, sigma );
    S=E(:,x,y);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    T=(DW_scheme1.sh'*DW_scheme1.sh+gamma*L);
    s=DW_scheme1.sh'*S;
    C=T\s;
    S1 = (HR_scheme1.sh*C)';
    S1(S1<0)=1e-4;
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    %[D4_coef D4]=estimate_ODF_with_SDP(E(:,x,y),E0,Gradients(:,1:3),order,delta);
    %[D6_coef D6]=estimate_ODF_with_SDP(E(:,x,y),E0,Gradients(:,1:3),6,delta);
    
    [D4_coef D4]=estimate_tensorODF(E(:,x,y),E0,Gradients(:,1:3),order,delta);
    [D6_coef D6]=estimate_tensorODF(E(:,x,y),E0,Gradients(:,1:3),6,delta);
    
%     if VOXEL.M==3
%         keyboard
%     end

%     if x==2 & y==10
%         keyboard
%     end
        
    P1=cp_als(tensor(D4),2,'init','nvecs','tol',1e-8); % compute directions with CP_ALS
    P2=cp_als(tensor(D6),3,'init','nvecs','tol',1e-8); % compute directions with CP_ALS
    % if P.lambda(2)<0.2 no_of_fibers=1 % for the testing data, compute the
    trueODF=VOXEL.ODF_true(ODF_XYZ(:,1),ODF_XYZ(:,2),ODF_XYZ(:,3));
    % compue estimated ODF
    G1=constructMatrixOfMonomials(ODF_XYZ,order);
    G2=constructMatrixOfMonomials(ODF_XYZ,6);
    
    myODF1=G1*D4_coef;
    myODF1=myODF1-min(myODF1);
    myODF1=myODF1/sum(myODF1); 
    
    myODF2=G2*D6_coef;
    myODF2=myODF2-min(myODF2);
    myODF2=myODF2/sum(myODF2);
    
    ODF_dev1(x,y)=norm(myODF1-trueODF)^2/norm(trueODF)^2;
    ODF_dev2(x,y)=norm(myODF2-trueODF)^2/norm(trueODF)^2;
    
    % number of fibers
    [dev1 mean_dev1] = calc_deviation(P1{1},VOXEL.R); %,VOXEL.M);
    [dev2 mean_dev2] = calc_deviation(P2{2},VOXEL.R); %,VOXEL.M);
    % mean_dev=calc_mean_deviation(P,VOXEL.R,VOXEL.M); 
    % mean_dev=1/2*(dev(1)+dev(2));
    mean_d1(x,y)=mean_dev1;
    mean_d2(x,y)=mean_dev2;
    
    %%%%%%%%%% low rank spherical deconvolution %%%%%%%%%
    
    %DT=Estimate_DTI(E(:,x,y),E0,Gradients(:,1:3),Gradients(1,4));
    %[Ev ev]=svd(DT);
    
    a=rand(3,3);
    Ev=Stabilized_Gram_Schmidt(a);
    
    [weights fiber_directions final_func]=deconvolve_two_fibers(E(:,x,y)', kernel, centroids, tri_areas, Gradients(:,1:3), Ev(:,1:2), options); %E(:,x,y)'
    %[weights fiber_directions final_func]=deconvolve_two_fibers(E1', kernel1, centroids, tri_areas, HR_scheme.vert, Ev(:,1:2), options); %E(:,x,y)'
    
    [dev3 mean_dev3] = calc_deviation(fiber_directions,VOXEL.R);
    mean_d3(x,y)=mean_dev3;
    %[fiber_direction iter]=low_rank_SD_one_fiber(E(:,x,y)', Ev(:,1), kernel, centroids, tri_areas, options.order, options.delta, options.tol);
    [weights fiber_directions_2 iter final_func] = deconvolve_three_fibers(E(:,x,y)', kernel, centroids, tri_areas, HR_scheme.vert, fiber_directions, options);
    [dev4 mean_dev4] = calc_deviation(fiber_directions_2,VOXEL.R);
    mean_d4(x,y)=mean_dev4;
    
    [ F_SH_srcsd, num_it ] = csdeconv (R_RH, E(:,x,y), DW_scheme, HR_scheme);
    %plot_SH (F_SH_srcsd, pscheme)
    Recon_ODF = SH2amp (F_SH_srcsd, pscheme);
    
    if VOXEL.M==1 | VOXEL.M==2
        P = find_max(Recon_ODF,pscheme.vert);
    elseif VOXEL.M==3
        P = find_max_3(Recon_ODF,pscheme.vert);
    end
    
    [dev5 mean_dev5] = calc_deviation(P.directions,VOXEL.R);   
    mean_d5(x,y)=mean_dev5;
    no_of_fibers(x,y)=VOXEL.M;
%     myODF3=0.5*(fiber_directions(:,1)'*ODF_XYZ').^16+0.5*(fiber_directions(:,2)'*ODF_XYZ').^16;
%     myODF3=myODF3-min(myODF3);
%     myODF3=myODF3/sum(myODF3);
%     ODF_dev3(x,y)=norm(myODF3'-trueODF)^2/norm(trueODF)^2;
    %[dev1 mean_dev1] = calc_deviation(fiber_directions,VOXEL.R(:,:,1:VOXEL.M),VOXEL.M);
    end
end

%%%%%%%%%%%%%%%% analyze results %%%%%%%%%%%%%%%%%%%%%%%%%
%[a1_LRD std1_LRD a2_LRD std2_LRD a3_LRD std3_LRD] = compute_average_deviation(mean_d4,no_of_fibers);
%[a1_CSD std1_CSD a2_CSD std2_CSD a3_CSD std3_CSD] = compute_average_deviation(mean_d5,no_of_fibers);

%save Structured_field_3D_results_SNR40_CSD_vs_LORPA_order12


