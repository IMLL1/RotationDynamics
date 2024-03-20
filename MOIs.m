%% ICs and propagate
clc; clear; close all;
rotm =  [1. 0 0; 0 0.866 -0.5;0 0.5 0.866];
points = [0   1  1  1;...
          0  -1  0  1;...
          0   1  0  1];
m = [1  1  1 1];
CoM = sum(m .* points,2)./sum(m);
points = points - CoM;
x = points(1,:); y = points(2, :); z = points(3, :);
Ix = sum(m.*(y.^2+z.^2)); Iy = sum(m.*(x.^2+z.^2)); Iz = sum(m.*(x.^2+y.^2));
Ixy = sum(m.*x.*y); Ixz = sum(m.*x.*z); Iyz = sum(m.*y.*z);
I = [Ix, Ixy, Ixz;Ixy, Iy, Iyz;Ixz, Iyz, Iz];
[V, D] = eigs(I);
PMOIs = diag(V'* I *V);
pointsNew = V'*points;
pointsNew = [pointsNew pointsNew(:,1)];
plot3(pointsNew(1,:),pointsNew(2,:),pointsNew(3,:)); axis equal; grid;