%% fem_postprocess.m
% Statistics and figures from fem_result.mat (written by fem_halbach_pde.m).
% Can be re-run on its own without re-solving.

here = fileparts(mfilename('fullpath'));
S = load(fullfile(here, 'fem_result.mat'));
fprintf('FEM mesh: hDSV=%.1f mm, hMag=%.1f mm, hAir=%.0f mm, mu_r=%.2f, Br=%.2f T\n', ...
    S.hDSV*1e3, S.hMag*1e3, S.hAir*1e3, S.muRec, S.Br);

Bmag = vecnorm(S.dsvB, 2, 2);
B0 = norm(S.B0vec);
P = S.dsvP;
fprintf('\n===== FEM RESULTS =====\n');
fprintf('B at centre       : [%+.4f %+.4f %+.4f] T, |B0| = %.2f mT\n', S.B0vec, B0*1e3);
fprintf('Mean |B| in target: %.2f mT  (1H: %.3f MHz)\n', mean(Bmag)*1e3, 42.577478*mean(Bmag));
fprintf('Raw |B| p-p in 40x40 cylinder : %.0f ppm (dominated by mesh noise)\n', ...
    (max(Bmag) - min(Bmag)) / mean(Bmag) * 1e6);
fprintf('Raw |B| RMS in 40x40 cylinder : %.0f ppm\n', std(Bmag) / mean(Bmag) * 1e6);

% Linear Nedelec elements give a piecewise-constant B, so raw point values
% are dominated by element-to-element jumps. The true field in the source-free
% bore is smooth, so fit |B| with a low-order 3D polynomial (least squares
% over many sample points) and report the fitted field's uniformity, plus
% the fit residual as the mesh-noise level.
Pn = P / (S.dsvDiam/2);
fitDeg = 4;
[ex, ey, ez] = ndgrid(0:fitDeg);
pw = [ex(:) ey(:) ez(:)];
pw = pw(sum(pw, 2) <= fitDeg, :);
Vm = ones(size(Pn, 1), size(pw, 1));
for j = 1:size(pw, 1)
    Vm(:, j) = Pn(:,1).^pw(j,1) .* Pn(:,2).^pw(j,2) .* Pn(:,3).^pw(j,3);
end
coef = Vm \ Bmag;
Bfit = Vm * coef;
noise = std(Bmag - Bfit) / mean(Bfit) * 1e6;
fprintf('\nSmoothed (degree-%d polynomial fit over %d points):\n', fitDeg, numel(Bmag));
fprintf('  fitted |B| at centre        : %.2f mT\n', coef(all(pw == 0, 2)) * 1e3);
fprintf('  fitted |B| p-p in 40x40 cyl : %.0f ppm\n', (max(Bfit) - min(Bfit)) / mean(Bfit) * 1e6);
fprintf('  fitted |B| RMS in 40x40 cyl : %.0f ppm\n', std(Bfit) / mean(Bfit) * 1e6);
fprintf('  mesh noise (fit residual RMS): %.0f ppm\n', noise);
for d = [20 30] * 1e-3
    k = hypot(P(:,1), P(:,2)) <= d/2 + 1e-12 & abs(P(:,3)) <= d/2 + 1e-12;
    fprintf('  fitted |B| p-p in %2.0f mm cyl : %.0f ppm\n', d*1e3, (max(Bfit(k)) - min(Bfit(k))) / mean(Bfit(k)) * 1e6);
end

% Figure 1: |B| maps in the central planes
s = S.planeS * 1e3; nS = numel(s);
f1 = figure('Color', 'w', 'Position', [50 50 1400 600]);
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; imagesc(s, s, reshape(vecnorm(S.planeXY, 2, 2), nS, nS)'); set(gca, 'YDir', 'normal');
axis equal tight; clim([0 1.5*B0]); colormap(turbo); cb = colorbar; cb.Label.String = '|B| [T]';
xlabel('x [mm]'); ylabel('y [mm]'); title('FEM |B|, z = 0');
nexttile; imagesc(s, s, reshape(vecnorm(S.planeXZ, 2, 2), nS, nS)'); set(gca, 'YDir', 'normal');
axis equal tight; clim([0 1.5*B0]); colormap(turbo); cb = colorbar; cb.Label.String = '|B| [T]';
xlabel('x [mm]'); ylabel('z [mm]'); title('FEM |B|, y = 0');
exportgraphics(f1, fullfile(here, 'fem_fig1_field_maps.png'), 'Resolution', 150);

% Figure 2: raw and smoothed profiles along the axes
t = S.lineT;
f2 = figure('Color', 'w', 'Position', [80 80 1300 500]);
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; hold on;
for i = 1:3, plot(t*1e3, vecnorm(S.lineB{i}, 2, 2)*1e3, 'LineWidth', 1); end
xline([-20 20], ':'); grid on; xlabel('position [mm]'); ylabel('|B| [mT]');
legend('along x', 'along y', 'along z', 'Location', 'best'); title('Raw FEM |B| along axes');
nexttile; hold on;
k = abs(t) <= S.dsvDiam/2;
axv = {[1 0 0], [0 1 0], [0 0 1]};
for i = 1:3
    q = (t(k) * axv{i}) / (S.dsvDiam/2);
    Vq = ones(nnz(k), size(pw, 1));
    for j = 1:size(pw, 1), Vq(:, j) = q(:,1).^pw(j,1) .* q(:,2).^pw(j,2) .* q(:,3).^pw(j,3); end
    plot(t(k)*1e3, (Vq*coef - coef(all(pw == 0, 2))) / coef(all(pw == 0, 2)) * 1e6, 'LineWidth', 1.5);
end
grid on; xlabel('position [mm]'); ylabel('deviation [ppm]');
legend('along x', 'along y', 'along z', 'Location', 'best');
title(sprintf('Smoothed FEM deviation (mesh noise ~%.0f ppm RMS)', noise));
exportgraphics(f2, fullfile(here, 'fem_fig2_profiles.png'), 'Resolution', 150);

% Figure 5: same layout as fig5_profiles.png from simulate_halbach_field.m,
% using the raw FEM field (lines were sampled over +/-25 mm only)
f5 = figure('Color', 'w', 'Position', [140 140 1100 450]);
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; hold on;
for i = 1:3, plot(t*1e3, vecnorm(S.lineB{i}, 2, 2)*1e3, 'LineWidth', 1.5); end
xline([-20 20], ':'); grid on; xlabel('position [mm]'); ylabel('|B| [mT]');
legend('along x', 'along y', 'along z', 'Location', 'best'); title('FEM |B| along the axes');
nexttile; hold on;
for i = 1:3, plot(t*1e3, (vecnorm(S.lineB{i}, 2, 2) - B0) / B0 * 1e6, 'LineWidth', 1.5); end
xline([-20 20], ':'); grid on; xlabel('position [mm]'); ylabel('deviation [ppm]');
legend('along x', 'along y', 'along z', 'Location', 'best');
title(sprintf('FEM deviation from B_0 = %.2f mT near centre', B0*1e3));
exportgraphics(f5, fullfile(here, 'fem_fig5_profiles.png'), 'Resolution', 150);
