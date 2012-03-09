function parResults (phantomName, bVal, SNR, denoise)

folderName = sprintf('Results/%s__B=%04d__SNR=%02d__DELTA=auto__DENOISE=%d', phantomName, bVal, SNR, denoise);
load( sprintf('TrainingData/%s__B=%04d__SNR=%02d__SIGNAL.mat', phantomName, bVal, SNR) );
[n1,n2,n3,nGrads] = size( E );

myESTIMATION = [];
myESTIMATION.FIELD = cell(  n1, n2, n3 );
myESTIMATION.ODF   = zeros( n1, n2, n3, 724 );

for x = 1:n1
for y = 1:n2
for z = 1:n3
    volFileName = fullfile(folderName, sprintf('v_%02d_%02d_%02d.mat', x,y,z));
    load (volFileName);
    myESTIMATION.FIELD{x,y,z} = result;    
    myESTIMATION.ODF(x,y,z,:) = getODF(result.R, result.lambda, result.f);
end
end
end

save( sprintf('Results/%s__B=%04d__SNR=%02d__DELTA=auto__DENOISE=%d.mat', phantomName, bVal, SNR, denoise), 'myESTIMATION' );

performance(phantomName, bVal, SNR, denoise);

end