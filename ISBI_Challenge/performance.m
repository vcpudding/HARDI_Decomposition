function performance (phantomName, bVal, SNR, denoise)

addpath HARDI_contest
load (sprintf('Results/%s__B=%04d__SNR=%02d__DELTA=auto__DENOISE=%d.mat', phantomName, bVal, SNR, denoise));
load( fullfile('TrainingData/Phantoms',phantomName) );
trueField = FIELD;
load( sprintf('TrainingData/Phantoms/%s_ODF.mat',phantomName) );
trueODF = ODF;

[n1,n2,n3] = size( trueField );

z = floor((n3+1)/2);
deviateBuf = zeros(n1, n2);
numDiffBuf = zeros(n1, n2);
odfDiffBuf = zeros(n1, n2);
numFiberBuf = zeros(n1, n2);

for x = 1:n1
for y = 1:n2
    nFibers = myESTIMATION.FIELD{x,y,z}.M;
    dirs = reshape(myESTIMATION.FIELD{x,y,z}.R(:, 3, :), [3, nFibers]);
    odf = myESTIMATION.ODF(x,y,z,:);
    
    deviateBuf(x,y) = mean(calc_deviation(dirs, trueField{x,y,z}.R));
    numDiffBuf(x,y) = nFibers-trueField{x,y,z}.M;
    odfDiffBuf(x,y) = norm(reshape(odf-trueODF(x,y,z,:), [724,1]))/norm(reshape(trueODF(x,y,z,:), [724,1]));
    numFiberBuf(x,y) = trueField{x,y,z}.M;
end
end

hFig = figure;
subplot(1,3,1); 
imagesc(deviateBuf); caxis([0, 20]); axis square; colorbar; 
title(sprintf('Direction deviation\n(%0.4f, %0.4f)\n(%0.4f, %0.4f)\n(%0.4f, %0.4f)', ...
    meanStat(deviateBuf, numFiberBuf, 1), stdStat(deviateBuf, numFiberBuf, 1), ...
    meanStat(deviateBuf, numFiberBuf, 2), stdStat(deviateBuf, numFiberBuf, 2), ...
    meanStat(deviateBuf, numFiberBuf, 3), stdStat(deviateBuf, numFiberBuf, 1)));
subplot(1,3,2); imagesc(numDiffBuf); caxis([-2, 2]); axis square; colorbar;
title(sprintf('Difference in # of fibers\n(%0.4f, %0.4f)\n(%0.4f, %0.4f)\n(%0.4f, %0.4f)', ...
    meanStat(numDiffBuf, numFiberBuf, 1), stdStat(numDiffBuf, numFiberBuf, 1), ...
    meanStat(numDiffBuf, numFiberBuf, 2), stdStat(numDiffBuf, numFiberBuf, 2), ...
    meanStat(numDiffBuf, numFiberBuf, 3), stdStat(numDiffBuf, numFiberBuf, 1)));
subplot(1,3,3); imagesc(odfDiffBuf); caxis([0, 0.5]); axis square; colorbar;
title(sprintf('Difference in ODF\n(%0.4f, %0.4f)\n(%0.4f, %0.4f)\n(%0.4f, %0.4f)', ...
    meanStat(odfDiffBuf, numFiberBuf, 1), stdStat(odfDiffBuf, numFiberBuf, 1), ...
    meanStat(odfDiffBuf, numFiberBuf, 2), stdStat(odfDiffBuf, numFiberBuf, 2), ...
    meanStat(odfDiffBuf, numFiberBuf, 3), stdStat(odfDiffBuf, numFiberBuf, 1)));
saveas(hFig, sprintf('../Data/figures/%s__B=%04d__SNR=%02d__DENOISE=%d.fig', phantomName, bVal, SNR, denoise));
%print(hFig, sprintf('../Data/figures/%s__B=%04d__SNR=%02d__DENOISE=%d.eps', phantomName, bVal, SNR, denoise), '-depsc');

close(hFig);

end

function y=meanStat(statBuf, numBuf, num)
numStatBuf = statBuf(numBuf==num);
y = mean(numStatBuf(:));
end

function y=stdStat(statBuf, numBuf, num)
numStatBuf = statBuf(numBuf==num);
y = std(numStatBuf(:));
end