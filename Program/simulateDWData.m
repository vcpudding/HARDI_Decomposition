function [S trueDirs] = simulateDWData (bVal, gradientDirections, fibDirs, weights, isAnisotropic)

if isAnisotropic
    L = diag([ 3e-4, 3e-4, 1.7e-3 ]);
%     L = diag([ 1.7e-3, 3e-4, 3e-4 ]);
else
    L = diag([7e-4, 7e-4, 7e-4]);
end

S = zeros(length(gradientDirections), 1);

azimuth = pi/3;
zenith  = pi/3;
M = [ cos(azimuth) -sin(azimuth) 0 ; sin(azimuth) cos(azimuth) 0 ; 0 0 1 ] * ...
    [ cos(zenith) 0 sin(zenith) ; 0 1 0 ; -sin(zenith) 0 cos(zenith) ];
initDir = M*[0 0 1]';
L = M*L*M';

for i=1:length(gradientDirections)
    for j=1:length(fibDirs)
        angle = fibDirs(j);
        R = [cos(angle) sin(angle) 0;-sin(angle) cos(angle) 0;0 0 1];
        D = R*L*R';
        S(i) = S(i) + weights(j)*exp(-bVal*(gradientDirections(i,:)*D*gradientDirections(i,:)'));
    end
end

trueDirs = zeros(3, length(fibDirs));
for j=1:length(fibDirs)
    angle = fibDirs(j);
    R = [cos(angle) sin(angle) 0;-sin(angle) cos(angle) 0;0 0 1];
    trueDirs(:,j) = R*initDir;
end

end