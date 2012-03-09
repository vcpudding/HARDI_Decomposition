function [dirs, weights, stat] = estimateFibers (S, gradientDirs, bVal, options, deconvOptions)

energyThres = options.energy_thres;
minNumValid = options.min_num_valid;
minNumTrial = options.min_num_trial;
maxNumTrial = options.max_num_trial;

nValid = 0;
dirBuf = [];
weightBuf = [];
energyBuf = [];
nFibers = 3;

for i=1:maxNumTrial
    if options.algorithm == 0
        [d, w, stat] = deconvolveFibersGD(S, gradientDirs, bVal, 3, deconvOptions);
    else
        [d, w, stat] = deconvolveFibersLM(S, gradientDirs, bVal, 3, deconvOptions);
    end
    
    if stat.convEnergy < energyThres
        d = reshape(d, [1,9]);
        n = length(w(w>options.fraction_thres));
        if n<nFibers
            dirBuf = d;
            weightBuf = w;
            energyBuf = stat.convEnergy;
            nFibers = n;
        else
            if n==nFibers
                dirBuf = [dirBuf; d];
                weightBuf = [weightBuf; w];
                energyBuf = [energyBuf; stat.convEnergy];
            end
        end
        nValid = nValid+1;
    end
    
    if (nValid>=minNumValid && i>=minNumTrial) || maxNumTrial-i<minNumValid-nValid
        break;
    end
end

if nValid<minNumValid
    dirs = [];
    weights = [];
    return;
else
    [minEnergy minIdx] = min(energyBuf);
    dirs = reshape(dirBuf(minIdx,1:3*nFibers), [3, nFibers]);
    weights = weightBuf(minIdx,1:nFibers);
%     dir1 = reshape(dirBuf(1,:), [3, 3]);
%     disp(dirBuf);
% 
%     for i=2:size(dirBuf,1)
%         dir2 = reshape(dirBuf(i,:), [3, 3]);
%         weight2 = weightBuf(i,:);
%         
%         for j=1:nFibers
%             d = dir1(:,j)'*dir2(:,j:nFibers);
%             [junk, idx] = max(abs(d));
%             idx = idx+j-1;
%             if d(idx-j+1)<0
%                 dir2(:,idx) = -dir2(:,idx);
%             end
%             tmp = dir2(:,idx);
%             dir2(:,idx) = dir2(:,j);
%             dir2(:,j) = tmp;
%             tmp = weight2(idx);
%             weight2(idx) = weight2(j);
%             weight2(j) = tmp;
%         end
%     
%         dirBuf(i,:) = reshape(dir2, [1, 9]);
%         weightBuf(i,:) = weight2;
%     end
%     disp(dirBuf);
%     
% %     dirs = mean(dirBuf(:, 1:3*nFibers), 1);
% %     weights = mean(weightBuf(:, 1:nFibers), 1);
%     dirs = zeros(3, nFibers);
%     weights = zeros(1, nFibers);
%     medIdx = floor((size(dirBuf,1)+1)/2);
%     for i=1:nFibers
%         [dirBuf sortIdx] = sortrows(dirBuf, 3*(i-1)+1);
%         dirs(:,i) = dirBuf(medIdx, 3*(i-1)+1:3*i)';
%         dirs(:,i) = dirs(:,i)/norm(dirs(:,i));
%         weights(i) = weightBuf(sortIdx(medIdx), i);
%     end
%     
%     weights = weights/sum(weights);
end

stat.convEnergy = min(energyBuf);

end