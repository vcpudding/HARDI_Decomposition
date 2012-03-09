function odf=getODF(R, lambda, f)
load( 'TrainingData/ODF_XYZ.mat' );

M = length(f);

invDBuf = zeros(size(R));
detDBuf = zeros(1, M);

for j=1:M
    D = R(:,:,j)'*diag(lambda(:,j))*R(:,:,j);
    invDBuf(:,:,j) = inv(D);
    detDBuf(j) = sqrt(det(D));
end

odf = zeros(724, 1);
for i=1:724
    dir = ODF_XYZ(i,:);
    for j=1:M
        odf(i) = odf(i)+f(j)*(dir*invDBuf(:,:,j)*dir')^(-3/2)/(4*pi*detDBuf(j));
    end
end
    
end