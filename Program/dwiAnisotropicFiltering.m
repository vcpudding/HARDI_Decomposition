function [filteredImg eBuf] = dwiAnisotropicFiltering (noisyImg, stdDeviation, conductance, priorWeight, tolerance, step)

%initialize
filteredImg = noisyImg;
sqStdDeviation = stdDeviation^2;
logPost = logPosterior(filteredImg, noisyImg, stdDeviation, conductance, priorWeight);
logPost = sum(logPost(:));
a=0.1;
b=0.5;
it = 1;
eBuf = [];
%gradient descent
while it<=200
    lastLogPost = logPost;
    
    x = filteredImg.*noisyImg/sqStdDeviation;
    B = -filteredImg/sqStdDeviation + besseli(1,x,1)./besseli(0,x,1).*noisyImg/sqStdDeviation;
    P = vectorAnisotropicDivergence(filteredImg, conductance);
    
    g = B+priorWeight*P;
    filteredImg = filteredImg + step*g;
    
    logPost = logPosterior(filteredImg, noisyImg, stdDeviation, conductance, priorWeight);
    logPost = sum(logPost(:));
    %disp(['posterior:', num2str(logPost)]);
    disp(['#', num2str(it), '  e:', num2str(logPost)]);
    eBuf = [eBuf, logPost];
    if logPost-lastLogPost <tolerance
        break;
    end;
    
    it = it+1;
    
end

%output
%figure, plot(eBuf);

end

function p = logPosterior(filteredImg, noisyImg, sqStdDeviation, conductance, priorWeight)

z = filteredImg.*noisyImg/sqStdDeviation;
p1 = log(noisyImg/sqStdDeviation)-(filteredImg.^2+noisyImg.^2)/2/sqStdDeviation+log(besseli(0,z,1))+abs(z);

imgGradSqNorm = vectorGradientSqNorm(filteredImg);
E = exp(-imgGradSqNorm/2/conductance);
p2 = priorWeight*conductance*E;

p = p1+p2;
%p=p2;

end

function g = vectorGradientSqNorm (vol)

vecdim = size(vol,4);
gx = zeros(size(vol));
gy = zeros(size(vol));
gz = zeros(size(vol));

for i=1:vecdim
    if size(vol,3)>1
        [g1 g2 g3] = gradient(vol(:,:,:,i));
        gx(:,:,:,i) = g1;
        gy(:,:,:,i) = g2;
        gz(:,:,:,i) = g3;
    else
        [g1 g2] = gradient(vol(:,:,1,i));
        gx(:,:,:,i) = g1;
        gy(:,:,:,i) = g2;
    end
end

g = gx.^2+gy.^2+gz.^2;

end


function div = vectorAnisotropicDivergence (vol, conductance)

[m n p vecdim] = size(vol);
div = zeros(size(vol));
for i=1:m
    for j=1:n
        for k=1:p
            gN = vol(max(1,i-1),j,k,:) - vol(i,j,k,:);
            gS = vol(min(m,i+1),j,k,:) - vol(i,j,k,:);
            gW = vol(i, max(1,j-1),k,:) - vol(i,j,k,:);
            gE = vol(i, min(n,j+1),k,:) - vol(i,j,k,:);
            gF = vol(i, j, max(1,k-1),:) - vol(i,j,k,:);
            gP = vol(i, j, min(p,k+1),:) - vol(i,j,k,:);
            
            div(i,j,k,:) = exp(-dot(gN,gN)/conductance)*gN + ...
                           exp(-dot(gS,gS)/conductance)*gS + ...
                           exp(-dot(gW,gW)/conductance)*gW + ...
                           exp(-dot(gE,gE)/conductance)*gE + ...
                           exp(-dot(gF,gF)/conductance)*gF + ...
                           exp(-dot(gP,gP)/conductance)*gP;
        end
    end
end
end


function div = vectorAnisotropicDivergenceScalar (vol, conductance)

[m n p vecdim] = size(vol);
div = zeros(size(vol));
for i=1:m
    for j=1:n
        for k=1:p
            gN = vol(max(1,i-1),j,k,:) - vol(i,j,k,:);
            gS = vol(min(m,i+1),j,k,:) - vol(i,j,k,:);
            gW = vol(i, max(1,j-1),k,:) - vol(i,j,k,:);
            gE = vol(i, min(n,j+1),k,:) - vol(i,j,k,:);
            gF = vol(i, j, max(1,k-1),:) - vol(i,j,k,:);
            gP = vol(i, j, min(p,k+1),:) - vol(i,j,k,:);
            
            div(i,j,k,:) = exp(-gN.^2/conductance).*gN + ...
                           exp(-gS.^2/conductance).*gS + ...
                           exp(-gW.^2/conductance).*gW + ...
                           exp(-gE.^2/conductance).*gE + ...
                           exp(-gF.^2/conductance).*gF + ...
                           exp(-gP.^2/conductance).*gP;
        end
    end
end
end