function [m1]=Stabilized_Gram_Schmidt(m2)

n=size(m2,2);
m1(:,1)=m2(:,1)/norm(m2(:,1)); % e1=u1/||u1||

for i=2:n,
    for j=1:(i-1),
      v=Proj(m2(:,j),m2(:,i));
      m2(:,i)=m2(:,i)-v; % stablized Gram Schmidt process
    end
    m1(:,i)=m2(:,i)/norm(m2(:,i)); % normalize the vectors
end


function  [v1]=Proj(v2,v3)

v1=dot(v2,v3)/norm(dot(v2,v2))*v2;



