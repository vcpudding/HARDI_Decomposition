function S=Simulate_DW_data(b,GradientDirections,orientation,angle,w1,w2)

% orientation is the first fiber orientation
% angle if the separation angle between the fibers

l=[1.7e-3,3e-4,3e-4]; % the diffusivities

R=[cos(angle) sin(angle) 0;-sin(angle) cos(angle) 0;0 0 1];
R1=[cos(orientation) sin(orientation) 0;-sin(orientation) cos(orientation) 0;0 0 1];
R2=R1*R;

L=diag(l);

D1=R1'*L*R1; % diffusion tensor 1
D2=R2'*L*R2; % diffusion tensor 2

% simulate data

for i=1:length(GradientDirections)
    S(i,1)=w1*exp(-b*GradientDirections(i,:)*D1*GradientDirections(i,:)')+w2*exp(-b*GradientDirections(i,:)*D2*GradientDirections(i,:)');
end



    

