function d = directionDeviation (dirs, trueDirs)

nDirs = size(dirs, 2);
nTrueDirs = size(trueDirs,2);
d = zeros(1, nDirs);
for i=1:nDirs
    d(i) = pi/2;
    for j=1:nTrueDirs
        dev = acos(dirs(:,i)'*trueDirs(:,j));
        if dev>pi/2
            dev = pi-dev;
        end
        d(i) = min(dev, d(i));
    end
end

d = d*180/pi;

end