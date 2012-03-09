
function [dev mean_dev] = calc_deviation(est, gt)

n=size(gt,3);
m=size(est,2);

for i=1:n
    for j=1:m
        dist_mat(i,j)=acos(abs(gt(:,3,i)'*est(:,j)))*180/pi;
    end
end

[x y]=assignmentoptimal(dist_mat);

mean_dev=y/n;

for k=1:n
 if x(k)~=0
     dev(k,1)=dist_mat(k,x(k));
 end
end

%        
% 
% deviation=[acos(abs(gt(:,3,1)'*est{1}(:,1)))*180/pi,acos(abs(gt(:,3,1)'*est{1}(:,2)))*180/pi,...
%                      acos(abs(gt(:,3,2)'*est{1}(:,1)))*180/pi,acos(abs(gt(:,3,2)'*est{1}(:,2)))*180/pi];
%                  
% [sorted idx]=sort(deviation);                 
%                  
% if idx(1)==1
%         dev(1)=acos(abs(gt(:,3,1)'*est{1}(:,1)))*180/pi;
%         dev(2)=acos(abs(gt(:,3,2)'*est{1}(:,2)))*180/pi;
%         
% elseif idx(1)==2
%         dev(1)=acos(abs(gt(:,3,1)'*est{1}(:,2)))*180/pi;
%         dev(2)=acos(abs(gt(:,3,2)'*est{1}(:,1)))*180/pi;
%         
% elseif idx(1)==3
%         dev(1)=acos(abs(gt(:,3,2)'*est{1}(:,1)))*180/pi;
%         dev(2)=acos(abs(gt(:,3,1)'*est{1}(:,2)))*180/pi;
%         
% else 
%         dev(1)=acos(abs(gt(:,3,2)'*est{1}(:,2)))*180/pi;
%         dev(2)=acos(abs(gt(:,3,1)'*est{1}(:,1)))*180/pi;
%         
% end
% 
% 
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% 
% 
% end

