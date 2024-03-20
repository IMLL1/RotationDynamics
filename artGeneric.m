%% ICs and propagate
clc; clear; close all;
dt = 1/100; tf = 15; % timestep and final time, both in sec
w0 = [.1;.2;10];  % initial angular velocity in the body frame

% MOIs
Ix = 100;
Iy = 500;
Iz = 300;

% Body shape, points are relative to CoM
points = [0     0  nan    -5     5;...
          0     0  nan     0     0;...
       -2.5  -0.5  nan  0.55  0.55];

% animation settings
tKeep = inf;         % time to keep traces
tracePts = [1, 4];  % which points to
traceColors = ['g', 'c'];
animSpd = 1.6;      % animation speed; to take lag into account

% Don't touch below
% Propagate
ops = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t,SV] = ode45(@(t,sv) dynamics (t,sv, Ix, Iy, Iz), 0:dt:tf, [0;0;0;w0], ops);
IC = "\omega_0^b=\langle"+w0(1)+", "+w0(2)+", "+w0(3)+"\rangle rad/s";

% calculate angular velocities and momenta and rotation matrices
wb = SV(:,4:6)'; % omega in body frame
rotms = eul2rotm(SV(:,1:3), "XYZ"); % rotation matrices for each frame
Hb = wb .* [Ix;Iy;Iz];
wi = reshape(wb, [3, 1, length(t)]);
wi = pagemtimes(rotms, wi); % omega inertial
wi = squeeze(wi);
Hi = reshape(Hb, [3, 1, length(t)]);
Hi = pagemtimes(rotms, Hi); % omega inertial
Hi = squeeze(Hi);

%% Euler Angles over time
eulSeqs = ["XZX", "XYX", "YXY","YZY" "ZYZ", "ZXZ", "XZY", "XYZ", "YXZ", "YZX", "ZYX", "ZXY"];
figure; tloEulerTime = tiledlayout(3,4, "TileSpacing","tight","Padding","tight");
title(tloEulerTime, "Euler Angles over Time - "+IC);
for idx = 1:length(eulSeqs)
    nexttile; euls = rotm2eul(rotms, eulSeqs(idx));
    plot(unwrap(euls), '-');
    actSeq = char(eulSeqs(idx));
    legend(["\phi: "+actSeq(1), "\theta: "+actSeq(2), "\psi: "+actSeq(3)], "Location","northwest", "Color","none", 'Box','off');
    xticklabels({});yticklabels({});
end
xlabel(tloEulerTime, "Time (t) [sec]");
ylabel(tloEulerTime, "Angle (\phi, \theta, \psi) [rad]");

%% Omega over time
figure; tloOmega = tiledlayout(1,2, "TileSpacing", "tight", "Padding", "tight");
title(tloOmega, "Angular Velocity - "+ IC);
nexttile; plot(t, wb); legend(["\omega_x^b", "\omega_y^b", "\omega_z^b"]);
title("Body Frame"); grid;
nexttile; plot(t, wi); legend(["\omega_x^i", "\omega_y^i", "\omega_z^i"]);
title("Inertial Frame"); grid on;
xlabel(tloOmega, "Time (t) [sec]");
ylabel(tloOmega, "Angular Velocity (\omega) [rad/s]");

%% omega, H in body and inertial (3d)
figure; maxw = max(vecnorm(wb)); maxH = max(vecnorm(Hb));
tloBodCone = tiledlayout(2,2, "TileSpacing","none","Padding","tight");
title(tloBodCone, "Angular Velocity and Momentum in 3D Space - "+IC);
ylabel(tloBodCone, "Angular Momentum - Angular Velocity");
xlabel(tloBodCone, "Body Frame - Inertial Frame");

nexttile; plot3(wb(1,:), wb(2,:), wb(3,:), '-b'); axis equal; grid on;
xlabel("x^b"); ylabel("y^b"); zlabel("z^b");
xlim([-maxw, maxw]*1.25); ylim([-maxw, maxw]*1.25); zlim([-maxw, maxw]*1.25);

nexttile; plot3(wi(1,:), wi(2,:), wi(3,:), '-b'); axis equal; grid on;
xlabel("x^i"); ylabel("y^i"); zlabel("z^i");
xlim([-maxw, maxw]*1.25); ylim([-maxw, maxw]*1.25); zlim([-maxw, maxw]*1.25);

nexttile; plot3(Hb(1,:), Hb(2,:), Hb(3,:), '-b'); axis equal; grid on;
xlabel("x^b"); ylabel("y^b"); zlabel("z^b");
xlim([-maxH, maxH]*1.25); ylim([-maxH, maxH]*1.25); zlim([-maxH, maxH]*1.25);

nexttile; plot3(Hi(1,:), Hi(2,:), Hi(3,:), 'ob'); axis equal; grid on;
xlabel("x^i"); ylabel("y^i"); zlabel("z^i");
xlim([-maxH, maxH]*1.25); ylim([-maxH, maxH]*1.25); zlim([-maxH, maxH]*1.25);

%% Attitude in 3d space
mappedX = squeeze(rotms(:,1,:));
mappedY = squeeze(rotms(:,2,:));
mappedZ = squeeze(rotms(:,3,:));

xX = mappedX(1,:); xY = mappedX(2,:); xZ = mappedX(3,:);
yX = mappedY(1,:); yY = mappedY(2,:); yZ = mappedY(3,:);
zX = mappedZ(1,:); zY = mappedZ(2,:); zZ = mappedZ(3,:);

figure; tloAttitude = tiledlayout(1,3, "TileSpacing","none","Padding","tight");
title(tloAttitude, "Unit Attitude In Inertial Space - " + IC);
subtitle(tloAttitude, IC+" for duration of "+tf+" sec"); l = [-1.25, 1.25];

nexttile; plot3(xX, xY, xZ, '-k'); grid on; axis equal; xlim(l); ylim(l); zlim(l);
xlabel("x^i [m]"); ylabel("y^i [m]"); zlabel("z^i [m]");
title("Unit Body x");

nexttile; plot3(yX, yY, yZ, '-k'); grid on; axis equal; xlim(l); ylim(l); zlim(l);
xlabel("x^i [m]"); ylabel("y^i [m]"); zlabel("z^i [m]");
title("Unit Body y");

nexttile; plot3(zX, zY, zZ, '-k'); grid on; axis equal; xlim(l); ylim(l); zlim(l);
xlabel("x^i [m]"); ylabel("y^i [m]"); zlabel("z^i [m]");
title("Unit Body z");

%% Euler Angles in 3d Space

eulSeqs = ["XZX", "XYX", "YXY","YZY" "ZYZ", "ZXZ", "XZY", "XYZ", "YXZ", "YZX", "ZYX", "ZXY"];
n = length(eulSeqs); r = floor(sqrt(n)); c = n/r;
figure; tloEul = tiledlayout(r, ceil(c), "TileSpacing","loose", "Padding","tight");
title(tloEul, "Euler Angles Over Time - "+IC);
for idx = 1:n
    nexttile; euls = rotm2eul(rotms, eulSeqs(idx));
    plot3(unwrap(euls(:,1)), unwrap(euls(:,2)), unwrap(euls(:,3)), '-k');
    grid on; axis equal; title(eulSeqs(idx));
    xlabel("\phi"); ylabel("\theta"); zlabel("\psi");
end
%% Tracked Points
pointsAllFrames = pagemtimes(rotms, repmat(points,1,1,size(rotms,3)));

figTracked = figure; plot3(nan,nan,nan); hold on; grid;
title("Tracked Points In Inertial Space - "+IC);

MS = max(vecnorm(points));
axis equal; xlim([-1, 1]*1.025*MS); ylim([-1, 1]*1.025*MS); zlim([-1, 1]*1.025*MS)
xticklabels({});yticklabels({});zticklabels({});

for idx=1:length(tracePts)
    X = squeeze(pointsAllFrames(1,tracePts(idx),:));
    Y = squeeze(pointsAllFrames(2,tracePts(idx),:));
    Z = squeeze(pointsAllFrames(3,tracePts(idx),:));
    plot3(X, Y, Z, "color", traceColors(idx), "LineWidth", 0.1);
end

%% Animation
nKeep = round(tKeep/dt);
pointsAllFrames = pagemtimes(rotms, repmat(points,1,1,size(rotms,3)));

animFig = figure('Color',[0.025 0.025 0.025]);  plot3(nan, nan, nan);
hold on; title("Animation In Inertial Space - "+IC, 'Color', 'w')%);
subtitle("Press any key to continue", 'Color', 'w')%);
for idx=1:length(tracePts)
    al(idx) = animatedline("Color",traceColors(idx), "LineWidth", 0.1, "MaximumNumPoints",nKeep);
end
obj = plot3(0,0,0, '-w', 'LineWidth',2); pens = plot3(0,0,0,'.', 'color', [1 0 0], 'MarkerSize',10);

MS = max(vecnorm(points));
plot3(MS*[-1.25, 1.25, nan, 0, 0, nan, 0, 0], MS*[0, 0, nan, -1.25, 1.25, nan, 0, 0], MS*[0, 0, nan, 0, 0, nan, -1.25, 1.25], '--w')
axis equal; xlim([-1, 1]*1.025*MS); ylim([-1, 1]*1.025*MS); zlim([-1, 1]*1.025*MS)
animAx = gca;
set(animAx,'Color',[0.025 0.025 0.025]); grid; set(animAx,'GridColor','w');
set(animAx,'XColor','w');set(animAx,'YColor','w');set(animAx,'ZColor','w'); xticklabels({});yticklabels({});zticklabels({});
pause; tic;
%animation = VideoWriter('animation'); animation.FrameRate = 1/dt; open(animation); % Comment to not save animation
for idx=1:length(t)
    points1frame = pointsAllFrames(:,:,idx);
    X = points1frame(1,:); Y = points1frame(2,:); Z = points1frame(3,:);
    set(obj, "XData", X, "YData", Y, "ZData", Z);
    set(pens, "XData", X(tracePts), "YData", Y(tracePts), "ZData", Z(tracePts))
    for trIdx=1:length(tracePts)
        addpoints(al(trIdx), X(tracePts(trIdx)), Y(tracePts(trIdx)), Z(tracePts(trIdx)));
    end
    animAx.Subtitle.String = "t = "+round(t(idx), 2)+" sec";
    %writeVideo(animation, getframe(gcf))
    pause(max([dt-toc, 0])/animSpd); tic;
end
%close(animation);

%% helper func
function dsv = dynamics(t, sv, Ix, Iy, Iz)
ph = sv(1); th = sv(2); ps = sv(3); % phi, theta, psi
wx = sv(4); wy = sv(5); wz = sv(6); % omega x, y, and z
dtheta = wx*sin(ps) + wy*cos(ps);
while(abs(cos(th)) <= 1e-8) % to avoid gimbal lock
    th = th + 1e-8 * dtheta;
end
dphi = (wx*cos(ps) - wy*sin(ps)) / cos(th);
dpsi = (-wx*cos(ps) + wy*sin(ps)) *sin(th) / cos(th) + wz;
dwx = (Iy-Iz)*wy*wz/Ix;
dwy = (Iz-Ix)*wx*wz/Iy;
dwz = (Ix-Iy)*wx*wy/Iz;
dsv = [dphi;dtheta;dpsi;dwx;dwy;dwz];
end