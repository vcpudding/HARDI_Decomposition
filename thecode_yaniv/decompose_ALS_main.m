
close all;
clear all;

%load sphere_triangulation_2562

%addpath /scratch/Documents/Postdoc/MATLAB/4thorder/Decomposition_project
addpath C:\Users\Yaniv\Documents\Postdoc\MATLAB\4thorder\Decomposition_project

UnitVectors;

g2=BuildSphere(3); % discrete integration directions: 2=162, 3=642, 4=2562 directions

GD=g(1:81,:); % gradient directions
b=3000; delta=500; order=4;
w1=0.5; w2=0.5; %.5;  %0.5;
t=pi/3; sig=0; %0.02; %0.04; %.05; % SNR=1/sig

epsilon=1e-6; %e-6;

r=2; % rank

fiber_orientation1=[0 1 0]; %/sqrt(2) 1/sqrt(2) 0];
orientation1=atan2(fiber_orientation1(2),fiber_orientation1(1));

R=[cos(t) -sin(t) 0;sin(t) cos(t) 0;0 0 1];
fiber_orientation2=fiber_orientation1*R';

Si=SimulateData(b,GD,t,orientation1,w1,w2);

S=Si;
y=randn(2,length(Si));
S = abs(Si+sig*(y(1,:)+sqrt(-1)*y(2,:)));

[SS idx] = sort(S);

%v_coef=rand(3,r);
v_coef(:,1)=[1/sqrt(2) 1/sqrt(2) 0]; %v_coef(:,1)/norm(v_coef(:,1));
v_coef(:,2)=[0 0 1]; %v_coef(:,2)/norm(v_coef(:,2));

%v_coef=[0 2 0 1 0 0]';

%v_coef(:,1)=g(idx(1),:); % initialize with the minimum of S
%v_coef(:,2)=g(idx(2),:);

v_coef_update=v_coef; %zeros(3,r);

dt=1;


% initialization step

% [dF M] = compute_gradient_ALS(S, v_coef_update, GD, order, delta, 1);
% Jcb=M*M';
% [ev eg]=eigs(Jcb);
% v_coef_update(:,1)=Jcb(:,1);
% v_coef_update(:,2)=Jcb(:,2);
tol=1e-6;
iter=1;
dFunc=1;
dNorm=1;

tic
while dNorm>tol

v=1;

v_coef_update_old=v_coef_update;

[dF M] = compute_gradient_ALS_fast(S, v_coef_update, g2, GD, order, delta, v); % area, cen

%[dF M H] = compute_gradient_ALS_2(S, v_coef_update, GD, order, delta, v); % area, cen

%F=dF.^2;
v_coef_update(:,1)=v_coef_update(:,1)-dt*inv(M'*M+epsilon*diag(diag(M'*M)))*M'*dF; % Newton step
%v_coef_update(:,1)=v_coef_update(:,1)-dt*inv(H)*M*dF'; % Newton step

%fnorm1=fnorm1*norm(v_coef_update(:,1));
%v_coef_update(:,1)=v_coef_update(:,1)/norm(v_coef_update(:,1));


v=2;
[dF M] = compute_gradient_ALS_fast(S, v_coef_update, g2, GD, order, delta, v); % area, cen

%[dF M H] = compute_gradient_ALS_2(S, v_coef_update, GD, order, delta, v); % area, cen

v_coef_update(:,2)=v_coef_update(:,2)-dt*inv(M'*M+epsilon*diag(diag(M'*M)))*M'*dF;
%v_coef_update(:,2)=v_coef_update(:,2)-dt*inv(H)*M*dF';

% normalize coefficients

fnorm1(iter)=norm(v_coef_update(:,1));
fnorm2(iter)=norm(v_coef_update(:,2));
%fnorm2=fnorm2*norm(v_coef_update(:,2));
%v_coef_update(:,1)=v_coef_update(:,1)/norm(v_coef_update(:,1));
%v_coef_update(:,2)=v_coef_update(:,2)/norm(v_coef_update(:,2));

F=0.5*dF.^2;
func(iter)=sum(F);
display(func(iter));

if iter>1
    FuncDiff=func(iter)-func(iter-1);
    dFunc=abs(FuncDiff)/func(iter); 
    dNorm=max(abs(fnorm1(iter)-fnorm1(iter-1)),abs(fnorm2(iter)-fnorm2(iter-1)));
    if FuncDiff<0 
        epsilon=epsilon/10;
        iter=iter+1;
    else
        epsilon=epsilon*10;
        v_coef_update=v_coef_update_old;
    end
else iter=iter+1;
end
    

if iter>400
    display('slow convergence');
    break
end

% else
%     v_coef_update=v_coef_update_temp;
%     iter=iter+1;
% end


%v_coef_update=v_coef_update-dt*M*F';
%v_coef=v_coef_update;

end

% refinment step
% for i=1:100
% 
% v=2;
% [dF M] = compute_gradient_ALS_fast(S, v_coef_update, g2, GD, order, delta, v); % area, cen
% 
% %v_coef_update(:,1)=v_coef_update(:,1)-dt*inv(M*M'+epsilon*diag(diag(M*M')))*M*dF';
% 
% v_coef_update(:,v)=v_coef_update(:,v)-dt*inv(M'*M+epsilon*diag(diag(M'*M)))*M'*dF;
% 
% %[dF M] = compute_gradient(S, v_coef_update(:,2), g2, GD, order, delta); % area, cen
% 
% %v_coef_update(:,2)=v_coef_update(:,2)-dt*inv(M*M'+epsilon*diag(diag(M*M')))*M*dF';
% end

toc
% fibers (rank-2)
v1=v_coef_update(:,1);
v2=v_coef_update(:,2);

fnorm1=norm(v1);
fnorm2=norm(v2);

display('ground-truth fiber directions');
[fiber_orientation1;fiber_orientation2]

display('estimated_fiber_directions');
fiber_direction1 = v1/fnorm1;
fiber_direction2 = v2/fnorm2;
[fiber_direction1';fiber_direction2']

display('estimated weights');
weight1 = fnorm1^order/(fnorm1^order+fnorm2^order);
weight2 = fnorm2^order/(fnorm1^order+fnorm2^order);
[weight1;weight2]

display('actual weights');
[w1;w2]

display(iter);
figure,plot(func);

