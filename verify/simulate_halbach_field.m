%% simulate_halbach_field.m
% Field simulation of the permanent-magnet array in
% BEST_bore100_cyl40x40_649mag_427ppm.xlsx
%
% Model: each magnet is a uniformly magnetized cuboid (exact analytic
% surface-charge solution, Furlani). Magnet-to-magnet demagnetization
% is ignored (mu_r ~ 1.05 for NdFeB, so this is a ~few-% effect on |B|).
%
% Rotation convention: Z-rot = 0 -> north (magnetization) points along +x;
% positive Z-rot rotates counter-clockwise about +z. X-rot / Y-rot are applied
% too (intrinsic Z*Y*X) but are all zero in this sheet.

clear; close all; clc;

%% ---------------- Magnet parameters (EDIT TO MATCH SPEC SHEET) -------------
Br        = 1.32;      % remanence [T]: N42, Br max 13,200 G per spec sheet
magSize   = [];        % [Lx Ly Lz] in mm, magnetization along local x.
                       % Empty -> cube of side Magnet_length from the sheet.
dsvDiam   = 40;        % target region: cylinder diameter [mm]
dsvLength = 40;        % target region: cylinder length along z [mm]
nGrid     = 41;        % grid points per axis inside the target region

%% ---------------- Load layout ---------------------------------------------
here = fileparts(mfilename('fullpath'));
T = readtable(fullfile(here, 'BEST_bore100_cyl40x40_649mag_427ppm.xlsx'), ...
              'VariableNamingRule', 'preserve');
pos = [T.('X-pos') T.('Y-pos') T.('Z-pos')];              % mm
rot = [T.('X-rot') T.('Y-rot') T.('Z-rot')];              % deg
if isempty(magSize), magSize = T.Magnet_length(1) * [1 1 1]; end
nMag = size(pos, 1);
R = zeros(3, 3, nMag);
for m = 1:nMag
    R(:,:,m) = Rz(rot(m,3)) * Ry(rot(m,2)) * Rx(rot(m,1));
end
fprintf('%d magnets, %.2f x %.2f x %.2f mm, Br = %.2f T\n', nMag, magSize, Br);

%% ---------------- Validate the cuboid model against a dipole sum ----------
validateCuboidModel(Br, magSize);

%% ---------------- Field in the target cylinder ----------------------------
g = linspace(-1, 1, nGrid);
[X, Y, Z] = ndgrid(g * dsvDiam/2, g * dsvDiam/2, g * dsvLength/2);
inDSV = hypot(X, Y) <= dsvDiam/2 + 1e-9;
P = [X(inDSV) Y(inDSV) Z(inDSV)];
B = arrayField(P, pos, R, magSize, Br);
Bmag = vecnorm(B, 2, 2);

B0vec = arrayField([0 0 0], pos, R, magSize, Br);
B0 = norm(B0vec);
[~, mainAx] = max(abs(B0vec));
axName = 'xyz';

ppmPP   = (max(Bmag) - min(Bmag)) / mean(Bmag) * 1e6;
ppmRMS  = std(Bmag) / mean(Bmag) * 1e6;
ppmMain = (max(B(:,mainAx)) - min(B(:,mainAx))) / abs(mean(B(:,mainAx))) * 1e6;

% Uniformity over a spherical DSV of the same diameter, for reference
inSph = vecnorm(P, 2, 2) <= dsvDiam/2 + 1e-9;
ppmSph = (max(Bmag(inSph)) - min(Bmag(inSph))) / mean(Bmag(inSph)) * 1e6;

fprintf('\n===== RESULTS =====\n');
fprintf('B at centre        : [%+.4f %+.4f %+.4f] T, |B0| = %.4f T (%.1f mT)\n', B0vec, B0, B0*1e3);
fprintf('Main field axis    : %s\n', axName(mainAx));
fprintf('Larmor freq (1H)   : %.3f MHz\n', 42.577478 * B0);
fprintf('Target cylinder    : %g mm dia x %g mm long, %d points\n', dsvDiam, dsvLength, size(P,1));
fprintf('  |B| min / mean / max : %.5f / %.5f / %.5f T\n', min(Bmag), mean(Bmag), max(Bmag));
fprintf('  |B| peak-to-peak  : %.0f ppm  (%.3f mT)\n', ppmPP, (max(Bmag)-min(Bmag))*1e3);
fprintf('  |B| RMS deviation : %.0f ppm\n', ppmRMS);
fprintf('  B%s peak-to-peak  : %.0f ppm\n', axName(mainAx), ppmMain);
fprintf('  Spectral spread   : %.1f kHz (p-p, 1H)\n', 42.577478e3 * (max(Bmag)-min(Bmag)));
fprintf('%g mm DSV sphere   : %.0f ppm p-p\n', dsvDiam, ppmSph);

% ppm over smaller concentric cylinders (dia = length = d)
fprintf('\nUniformity vs region size (cylinder, dia = length):\n');
for d = [10 20 30 40]
    k = hypot(P(:,1), P(:,2)) <= d/2 + 1e-9 & abs(P(:,3)) <= d/2 + 1e-9;
    fprintf('  %2d mm : %6.0f ppm p-p, %6.0f ppm RMS\n', d, ...
        (max(Bmag(k))-min(Bmag(k)))/mean(Bmag(k))*1e6, std(Bmag(k))/mean(Bmag(k))*1e6);
end

%% ---------------- Figure 1: 3D layout + field slices ----------------------
nS = 61;  L = 95;  % slice extent [mm]
s = linspace(-L, L, nS);
[Sx, Sy] = ndgrid(s, s);
Bxy = reshape(vecnorm(arrayField([Sx(:) Sy(:) 0*Sx(:)], pos, R, magSize, Br), 2, 2), nS, nS);
Bxz = reshape(vecnorm(arrayField([Sx(:) 0*Sx(:) Sy(:)], pos, R, magSize, Br), 2, 2), nS, nS);
Byz = reshape(vecnorm(arrayField([0*Sx(:) Sx(:) Sy(:)], pos, R, magSize, Br), 2, 2), nS, nS);

f1 = figure('Name', '3D layout and |B|', 'Color', 'w', 'Position', [50 50 1100 850]);
hold on;
drawMagnets(pos, R, magSize);
surf(Sx, Sy, 0*Sx, Bxy, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
surf(Sx, 0*Sx, Sy, Bxz, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
surf(0*Sx, Sx, Sy, Byz, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
drawCylinder(dsvDiam/2, dsvLength);
colormap(f1, turbo); clim([0 1.5*B0]);
cb = colorbar; cb.Label.String = '|B| [T]';
axis equal; grid on; view(35, 22);
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
title(sprintf('Halbach array (%d magnets): |B|, centre = %.1f mT', nMag, B0*1e3));
lighting gouraud; camlight;
exportgraphics(f1, fullfile(here, 'fig1_layout_3d.png'), 'Resolution', 150);

%% ---------------- Figure 2: in-bore field lines ---------------------------
nQ = 15; q = linspace(-45, 45, nQ);
[Qx, Qy] = ndgrid(q, q);
keep = hypot(Qx, Qy) < 48;
Pq = [Qx(keep) Qy(keep) 0*Qx(keep)];
Bq = arrayField(Pq, pos, R, magSize, Br);

f2 = figure('Name', 'Field direction', 'Color', 'w', 'Position', [80 80 900 800]);
hold on;
imagesc(s, s, Bxy'); set(gca, 'YDir', 'normal');
quiver(Pq(:,1), Pq(:,2), Bq(:,1), Bq(:,2), 0.6, 'k', 'LineWidth', 1);
z0 = abs(pos(:,3)) < 1e-6;
mdir = squeeze(R(:,1,z0))';
quiver(pos(z0,1), pos(z0,2), mdir(:,1), mdir(:,2), 0.6, 'r', 'LineWidth', 1.2, 'MaxHeadSize', 0.8);
rectangle('Position', [-50 -50 100 100], 'Curvature', 1, 'LineStyle', '--');
colormap(f2, turbo); clim([0 1.5*B0]);
cb = colorbar; cb.Label.String = '|B| [T]';
axis equal tight; xlabel('x [mm]'); ylabel('y [mm]');
title('z = 0 plane: |B| (colour), field direction (black), magnet north (red)');
exportgraphics(f2, fullfile(here, 'fig2_field_direction_z0.png'), 'Resolution', 150);

%% ---------------- Figure 3: homogeneity maps (ppm) ------------------------
nH = 81; h = linspace(-dsvDiam/2, dsvDiam/2, nH);
[H1, H2] = ndgrid(h, h);
planes = {'xy (z=0)', [H1(:) H2(:) 0*H1(:)], 'x', 'y';
          'xz (y=0)', [H1(:) 0*H1(:) H2(:)], 'x', 'z';
          'yz (x=0)', [0*H1(:) H1(:) H2(:)], 'y', 'z'};
f3 = figure('Name', 'Homogeneity', 'Color', 'w', 'Position', [100 100 1500 480]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact');
for i = 1:3
    nexttile;
    ppm = (vecnorm(arrayField(planes{i,2}, pos, R, magSize, Br), 2, 2) - B0) / B0 * 1e6;
    ppm = reshape(ppm, nH, nH);
    if i == 1, ppm(hypot(H1, H2) > dsvDiam/2) = NaN; end
    contourf(h, h, ppm', 30, 'LineColor', 'none'); hold on;
    contour(h, h, ppm', [-400 -200 -100 -50 0 50 100 200 400], 'k', 'ShowText', 'on');
    axis equal tight; colorbar;
    xlabel([planes{i,3} ' [mm]']); ylabel([planes{i,4} ' [mm]']);
    title(['(|B| - B_0)/B_0 [ppm], ' planes{i,1}]);
end
colormap(f3, redblue(256));
title(tl, sprintf('Field deviation in %g mm target region, B_0 = %.2f mT', dsvDiam, B0*1e3));
exportgraphics(f3, fullfile(here, 'fig3_homogeneity_ppm.png'), 'Resolution', 150);

%% ---------------- Figure 4: 3D iso-surfaces of deviation ------------------
Dppm = nan(size(X)); Dppm(inDSV) = (Bmag - B0) / B0 * 1e6;
f4 = figure('Name', 'Deviation iso-surfaces', 'Color', 'w', 'Position', [120 120 900 800]);
hold on;
lv = max(abs(Dppm(:))) * [-0.5 0.5];
cols = [0.2 0.3 1; 1 0.25 0.2];
for i = 1:2
    p = patch(isosurface(permute(X,[2 1 3]), permute(Y,[2 1 3]), permute(Z,[2 1 3]), ...
                         permute(Dppm,[2 1 3]), lv(i)));
    set(p, 'FaceColor', cols(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.55);
end
hs = slice(permute(X,[2 1 3]), permute(Y,[2 1 3]), permute(Z,[2 1 3]), permute(Dppm,[2 1 3]), 0, 0, 0);
set(hs, 'EdgeColor', 'none', 'FaceColor', 'interp'); colormap(f4, redblue(256)); clim(max(abs(Dppm(:))) * [-1 1]);
cb = colorbar; cb.Label.String = 'deviation [ppm]';
drawCylinder(dsvDiam/2, dsvLength);
axis equal; grid on; view(35, 22); camlight; lighting gouraud;
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
title(sprintf('Deviation from B_0 in target cylinder; iso-surfaces at %+.0f / %+.0f ppm', lv));
exportgraphics(f4, fullfile(here, 'fig4_deviation_3d.png'), 'Resolution', 150);

%% ---------------- Figure 5: line profiles ---------------------------------
t = linspace(-50, 50, 401)';
f5 = figure('Name', 'Profiles', 'Color', 'w', 'Position', [140 140 1100 450]);
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; hold on;
lineP = {[t 0*t 0*t], [0*t t 0*t], [0*t 0*t t]};
for i = 1:3
    plot(t, vecnorm(arrayField(lineP{i}, pos, R, magSize, Br), 2, 2) * 1e3, 'LineWidth', 1.5);
end
xline([-20 20], ':'); grid on; xlabel('position [mm]'); ylabel('|B| [mT]');
legend('along x', 'along y', 'along z', 'Location', 'best'); title('|B| along the axes');
nexttile; hold on;
k = abs(t) <= 25;
for i = 1:3
    b = vecnorm(arrayField(lineP{i}(k,:), pos, R, magSize, Br), 2, 2);
    plot(t(k), (b - B0) / B0 * 1e6, 'LineWidth', 1.5);
end
xline([-20 20], ':'); grid on; xlabel('position [mm]'); ylabel('deviation [ppm]');
legend('along x', 'along y', 'along z', 'Location', 'best'); title('Deviation from B_0 near centre');
exportgraphics(f5, fullfile(here, 'fig5_profiles.png'), 'Resolution', 150);

fprintf('\nFigures written to %s\n', here);

%% ======================= Local functions ==================================
function M = Rx(d), M = [1 0 0; 0 cosd(d) -sind(d); 0 sind(d) cosd(d)]; end
function M = Ry(d), M = [cosd(d) 0 sind(d); 0 1 0; -sind(d) 0 cosd(d)]; end
function M = Rz(d), M = [cosd(d) -sind(d) 0; sind(d) cosd(d) 0; 0 0 1]; end

function B = arrayField(P, pos, R, dims, Br)
% Total B [T] at points P (N x 3, mm) from all magnets.
B = zeros(size(P));
for m = 1:size(pos, 1)
    Rm = R(:,:,m);
    Pl = (P - pos(m,:)) * Rm;              % world -> local:  Rm' * p
    Bl = cuboidFieldX(Pl, dims/2, Br);     % local field, M along local +x
    B = B + Bl * Rm';                      % local -> world: Rm * b
end
end

function B = cuboidFieldX(P, hd, Br)
% Exact field of a cuboid with half-dims hd = [a b c], magnetized along +x.
% Uses the z-magnetized solution with cyclic axis relabelling
% (x,y,z) -> (u,v,w) = (y,z,x), so magnetization lies along w.
Bw = cuboidFieldZ(P(:,[2 3 1]), hd([2 3 1]), Br);
B = Bw(:, [3 1 2]);
end

function B = cuboidFieldZ(P, hd, Br)
% Analytic field of a cuboid magnetized along +z, centred at 0 (surface
% charge +M on z=+c, -M on z=-c; Furlani). Returns mu0*(H+M), i.e. B,
% both outside and inside the magnet.
x = P(:,1); y = P(:,2); z = P(:,3);
a = hd(1); b = hd(2); c = hd(3);
xs = [-a a]; ys = [-b b]; zs = [-c c];
Hx = 0; Hy = 0; Hz = 0;
for i = 1:2, for j = 1:2, for k = 1:2
    dx = x - xs(i); dy = y - ys(j); dz = z - zs(k);
    r = sqrt(dx.^2 + dy.^2 + dz.^2);
    sgn = (-1)^(i + j + k);
    Hx = Hx - sgn * logDPlusR(dy, r, dx.^2 + dz.^2);
    Hy = Hy - sgn * logDPlusR(dx, r, dy.^2 + dz.^2);
    Hz = Hz + sgn * atan(dx .* dy ./ (dz .* r));
end, end, end
B = Br / (4*pi) * [Hx Hy Hz];
inside = abs(x) < a & abs(y) < b & abs(z) < c;
B(inside, 3) = B(inside, 3) + Br;
end

function v = logDPlusR(d, r, rho2)
% log(d + r) evaluated without cancellation when d ~ -r.
v = zeros(size(d));
p = d >= 0;
v(p) = log(d(p) + r(p));
v(~p) = log(max(rho2(~p), 1e-30) ./ (r(~p) - d(~p)));
end

function validateCuboidModel(Br, dims)
% Compare analytic cuboid to a fine sum of point dipoles (outside the magnet),
% and check B inside the magnet.
n = 24;
g = ((1:n) - 0.5) / n - 0.5;
[gx, gy, gz] = ndgrid(g*dims(1), g*dims(2), g*dims(3));
src = [gx(:) gy(:) gz(:)];
mdip = Br / (4*pi*1e-7) * prod(dims) * 1e-9 / numel(gx);   % A m^2 each, along x
tests = [1.2*dims(1) 0 0; 0 1.3*dims(2) 0; 0.9*dims(1) 0.8*dims(2) 0.7*dims(3); ...
         -2*dims(1) 1*dims(2) -1.5*dims(3); 0 0 3*dims(3)];
Ba = cuboidFieldX(tests, dims/2, Br);
Bd = zeros(size(tests));
for t = 1:size(tests, 1)
    rv = (tests(t,:) - src) * 1e-3;  rn = vecnorm(rv, 2, 2);
    Bd(t,:) = sum(1e-7 * mdip * (3 * rv(:,1) .* rv ./ rn.^5 - [1 0 0] ./ rn.^3), 1);
end
err = max(vecnorm(Ba - Bd, 2, 2) ./ vecnorm(Bd, 2, 2));
fprintf('Cuboid model vs dipole-sum check: max rel. error = %.2e\n', err);
if err > 1e-2
    error('Cuboid field model failed validation (rel err %.3g).', err);
end
end

function drawMagnets(pos, R, dims)
% Draw each magnet as a box coloured by magnetization direction angle.
[cx, cy, cz] = ndgrid([-1 1], [-1 1], [-1 1]);
corners = [cx(:) cy(:) cz(:)] .* dims / 2;
faces = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
V = zeros(8 * size(pos,1), 3); F = zeros(6 * size(pos,1), 4); C = zeros(8 * size(pos,1), 3);
for m = 1:size(pos, 1)
    V(8*(m-1)+(1:8), :) = corners * R(:,:,m)' + pos(m,:);
    F(6*(m-1)+(1:6), :) = faces + 8*(m-1);
    ang = atan2(R(2,1,m), R(1,1,m));
    C(8*(m-1)+(1:8), :) = repmat(hsv2rgb([mod(ang, 2*pi)/(2*pi) 0.7 0.9]), 8, 1);
end
patch('Vertices', V, 'Faces', F, 'FaceVertexCData', C, 'FaceColor', 'flat', ...
      'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.3);
mdir = squeeze(R(:,1,:))';
quiver3(pos(:,1), pos(:,2), pos(:,3), mdir(:,1), mdir(:,2), mdir(:,3), 0.5, 'k');
end

function drawCylinder(r, len)
th = linspace(0, 2*pi, 60);
for zz = [-len/2 len/2]
    plot3(r*cos(th), r*sin(th), zz + 0*th, 'k-', 'LineWidth', 1.5);
end
for a = 0:pi/2:3*pi/2
    plot3(r*cos(a)*[1 1], r*sin(a)*[1 1], [-len len]/2, 'k-', 'LineWidth', 1.5);
end
end

function cmap = redblue(n)
x = linspace(-1, 1, n)';
cmap = [min(1, 1 + x) .* (x < 0) + (x >= 0), 1 - abs(x), min(1, 1 - x) .* (x > 0) + (x <= 0)];
cmap = max(0, min(1, cmap));
end
