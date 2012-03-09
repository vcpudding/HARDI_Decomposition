function [fibDirs fibWeights stat] = deconvolveFibersGD (S, gradientDirs, bVal, nFibers, options)

order=options.order; % order=12 works very well
delta=options.delta; % estimated from the data
tol=options.tol; % set to 1e-6
innerTol = options.innertol;
maxIt=options.maxiter; % maximal number of iterations, default is 200
maxInnerIt = options.maxinneriter;
step = options.step;

if options.accurate_integration
    load sphere_integration_2562
else 
    load sphere_integration_642
end
triangleVecs = centroids'; % 3 x nVecs
triangleAreas = tri_areas; % 1 x nVecs


%%%%%%%%%%%%%%test delta%%%%%%%%%%%%%%%%%%%
% deltaBuf = 10:0.2:50;
% eBuf = zeros(size(deltaBuf));
% dBuf = zeros(size(deltaBuf));
% for i=1:length(deltaBuf)
%     delta = deltaBuf(i);
%     kernel = exp(-delta*(gradientDirs'*triangleVecs).^2);
%     fibDirs = [1 0 0; 0 1 0]';
%     energy = energyFunc(fibDirs, kernel, order, triangleVecs, triangleAreas, S);
%     eBuf(i) = energy;
%     dDir = dirGradient(fibDirs, kernel, order, triangleVecs, triangleAreas, S, 1);
%     dBuf(i) = norm(dDir(:));
% end
% figure;
% subplot(1,2,1); plot(deltaBuf, eBuf); title('Signal difference');
% subplot(1,2,2); plot(deltaBuf, dBuf); title('Derivative norm');
% fibWeights = ones(1,nFibers);
% stat = {};
% return;
%%%%%%%%%%%end of test delta%%%%%%%%%%%%%%%

kernel = exp(-delta*(gradientDirs'*triangleVecs).^2); % nGrads x nVecs

%initialization
switch options.init
    case 1
        % DTI
        DT=Estimate_DTI(S,1,gradientDirs',bVal);
        [Ev ev]=svd(DT);
        fibDirs = Ev(:,1:nFibers);
    case 2
        fibDirs = eye(3);
        fibDirs = fibDirs(:, 1:nFibers);
    otherwise
        fibDirs = rand(3, 3);
        fibDirs = Stabilized_Gram_Schmidt(fibDirs); % orthonormalize
        fibDirs = fibDirs(:, 1:nFibers);
end

energy = energyFunc(fibDirs, kernel, order, triangleVecs, triangleAreas, S);
disp(['initial energy: ', num2str(energy)]);
it = 0;
eBuf = zeros(maxIt,1);

a = 0.001;
b = 0.5;

stat.initEnergy = energy;
stat.initDirs = fibDirs;
tic;

while it<maxIt

    it = it+1;
    
    %update
    for i=1:nFibers
        
        if ~options.innerconvergence
            dDir = dirGradient(fibDirs, kernel, order, triangleVecs, triangleAreas, S, i);
            if options.linesearch
%                 disp('line search');
                dDir = dDir/norm(dDir(:,i));
                step = 1;
                while  step>1e-4 && energyFunc(fibDirs-step*dDir, kernel, order, triangleVecs, triangleAreas, S) > energyFunc(fibDirs, kernel, order, triangleVecs, triangleAreas, S)-a*step
                    step = b*step;            
                end
%                 disp(step);
            end
            fibDirs = fibDirs - step*dDir;
        else
            innerIt = 0;
            while innerIt<maxInnerIt
                innerIt = innerIt+1;
                dDir = dirGradient(fibDirs, kernel, order, triangleVecs, triangleAreas, S, i);
                if options.linesearch
                    step = 1;
                    while  step>1e-3 && energyFunc(fibDirs-step*dDir, kernel, order, triangleVecs, triangleAreas, S) > energyFunc(fibDirs, kernel, order, triangleVecs, triangleAreas, S)-a*step*norm(dDir(:,i))
                        step = b*step;            
                    end
                end
                fibDirs = fibDirs-step*dDir;
                e = norm(dDir(:,i));
                if  norm(e)<innerTol || isnan(e)
                    break;
                end
            end
        end
    end
    
%     disp(fibDirs);
    
    %check for termination
    lastEnergy = energy;
    energy = energyFunc(fibDirs, kernel, order, triangleVecs, triangleAreas, S);
    eBuf(it) = energy;
    if abs(energy-lastEnergy)<tol || isnan(energy)
        break;
    end
end

disp(['convergence energy: ', num2str(energy)]);

%finish up
fibNorms = sqrt(sum(fibDirs.^2));
tmpDirs = fibDirs./fibNorms(ones(3,1), :);
fibWeights = fibNorms.^order/sum(fibNorms.^order);
[fibWeights sortIdx] = sort(fibWeights, 'descend');
for i=1:nFibers
    fibDirs(:,i) = tmpDirs(:, sortIdx(i));
end

stat.convEnergy = energy;
stat.nIts = it;
stat.time = toc;
stat.eCurve = eBuf(1:it);

end

function e = simSigErr (fibDirs, kernel, order, triangleVecs, triangleAreas, S)

nFibers = size(fibDirs,2);
triangleFODs = ones(1, nFibers)*((fibDirs'*triangleVecs).^order); % 1 x nVecs
simS = kernel*(triangleFODs'.*triangleAreas'); % nGrads x 1
e = S-simS;

end

function e = energyFunc (fibDirs, kernel, order, triangleVecs, triangleAreas, S)

sigErr = simSigErr(fibDirs, kernel, order, triangleVecs, triangleAreas, S);
e = sum(sigErr(:).^2);

end

function dDir = dirGradient (fibDirs, kernel, order, triangleVecs, triangleAreas, S, updateFibIdx)

dDir = zeros(size(fibDirs));
sigErr = simSigErr(fibDirs, kernel, order, triangleVecs, triangleAreas, S);
fibDir = fibDirs(:, updateFibIdx);
dFODs = (fibDir'*triangleVecs).^(order-1); % 1 x nVecs
dFODxArea = dFODs.*triangleAreas; %1 x nVecs
c = -order*(dFODxArea(ones(3,1),:).*triangleVecs)*kernel'; % 3 x nGrads
dDir(:, updateFibIdx) = c*sigErr;

end
