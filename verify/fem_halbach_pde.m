%% fem_halbach_pde.m
% 3D finite-element magnetostatic simulation of the magnet array in
% BEST_bore100_cyl40x40_649mag_427ppm.xlsx using MATLAB's
% Partial Differential Equation Toolbox (femodel, AnalysisType="magnetostatic").
%
% Each magnet is a cell with a permanent magnetization load and recoil
% permeability mu_r, so magnet-magnet interaction (demagnetization) is
% included. Toolbox constitutive law:  B = mu0*mu_r*H + mu0*M,  M = Br/mu0.
%
% Rotation convention: Z-rot = 0 -> north along +x, CCW about +z.
% All geometry in SI units (m).
%
% Set meshOnly = true to just build/mesh and report the problem size.

if ~exist('meshOnly', 'var'), meshOnly = false; end
if ~exist('hDSV',   'var'), hDSV   = 2.0e-3; end   % mesh size on target-region surface [m]
if ~exist('hMag',   'var'), hMag   = 3.0e-3; end   % mesh size on magnet faces [m]
if ~exist('hAir',   'var'), hAir   = 30e-3;  end   % max mesh size far away [m]
if ~exist('hGrad',  'var'), hGrad  = 1.5;    end   % mesh growth rate

%% ---------------- Parameters ----------------------------------------------
Br        = 1.32;          % N42, Br max 13,200 G (spec sheet)
muRec     = 1.05;          % NdFeB recoil permeability
L         = 6.35e-3;       % cube edge [m]
dsvDiam   = 40e-3;         % target cylinder diameter [m]
dsvLength = 40e-3;         % target cylinder length [m]
rAir      = 0.20;          % air domain radius [m]
hAirLen   = 0.44;          % air domain length [m]
mu0       = 4*pi*1e-7;

here = fileparts(mfilename('fullpath'));
T = readtable(fullfile(here, 'BEST_bore100_cyl40x40_649mag_427ppm.xlsx'), ...
              'VariableNamingRule', 'preserve');
pos = [T.('X-pos') T.('Y-pos') T.('Z-pos')] * 1e-3;
th  = T.('Z-rot');
assert(all(T.('X-rot') == 0 & T.('Y-rot') == 0), 'Only Z rotations are handled.');
nMag = height(T);

%% ---------------- Geometry ------------------------------------------------
tStart = tic;
g = translate(fegeometry(multicylinder(rAir, hAirLen)), [0 0 -hAirLen/2]);
% Slightly oversized refinement cell around the target region
dsvCell = translate(fegeometry(multicylinder(dsvDiam/2 + 2e-3, dsvLength + 4e-3)), ...
                    [0 0 -(dsvLength + 4e-3)/2]);
g = addCell(g, dsvCell);
cube = translate(fegeometry(multicuboid(L, L, L)), [0 0 -L/2]);
for k = 1:nMag
    g = addCell(g, translate(rotate(cube, th(k), [0 0 0], [0 0 1]), pos(k,:)));
end
% addCell appends cells in insertion order: 1 = outer air, 2 = refinement
% cell, 3.. = magnets in spreadsheet order. Checked against the mesh below.
dsvCellId = 2;
magCell = (3:g.NumCells)';
assert(numel(magCell) == nMag);
fprintf('Geometry: %d cells (%.1f s)\n', g.NumCells, toc(tStart));

%% ---------------- Model ---------------------------------------------------
model = femodel(AnalysisType="magnetostatic", Geometry=g);
model.VacuumPermeability = mu0;
model.MaterialProperties = materialProperties(RelativePermeability=1);
model.MaterialProperties(magCell) = materialProperties(RelativePermeability=muRec);
for k = 1:nMag
    model.CellLoad(magCell(k)) = cellLoad(Magnetization=Br/mu0 * [cosd(th(k)); sind(th(k)); 0]);
end
innerFaces = unique(cell2mat(arrayfun(@(c) cellFaces(g, c), 2:g.NumCells, 'UniformOutput', false)));
outerFaces = setdiff(cellFaces(g, 1), innerFaces);
model.FaceBC(outerFaces) = faceBC(MagneticPotential=[0;0;0]);

magFaces = unique(cell2mat(arrayfun(@(c) cellFaces(g, c), magCell', 'UniformOutput', false)));
tMesh = tic;
model = generateMesh(model, GeometricOrder="linear", Hmax=hAir, Hgrad=hGrad, ...
    Hface={cellFaces(g, dsvCellId), hDSV, magFaces, hMag});
msh = model.Geometry.Mesh;
nEdgesEst = round(1.2 * size(msh.Elements, 2));   % ~edges per tet for Delaunay meshes
fprintf('Mesh: hDSV=%.1f mm, hMag=%.1f mm, hAir=%.0f mm -> %d nodes, %d tets (~%d edge DOFs) (%.1f s)\n', ...
    hDSV*1e3, hMag*1e3, hAir*1e3, size(msh.Nodes, 2), size(msh.Elements, 2), nEdgesEst, toc(tMesh));
% Verify every magnet cell sits where the spreadsheet says it should
chk = [1:25:nMag nMag];
for k = chk
    e = findElements(msh, 'region', 'Cell', magCell(k));
    cen = mean(msh.Nodes(:, unique(msh.Elements(:, e))), 2)';
    assert(norm(cen - pos(k,:)) < L/4, 'Cell %d is not magnet %d', magCell(k), k);
end
if meshOnly, return; end

%% ---------------- Solve ---------------------------------------------------
tSolve = tic;
R = solve(model);
fprintf('Solve: %.1f s\n', toc(tSolve));

%% ---------------- Post-processing -----------------------------------------
% The toolbox results object does not survive save/load (its interpolants are
% lost), so sample everything needed for plots now and save plain arrays.
dsvN = 41;   gv = linspace(-1, 1, dsvN);
[X, Y, Z] = ndgrid(gv * dsvDiam/2, gv * dsvDiam/2, gv * dsvLength/2);
in = hypot(X, Y) <= dsvDiam/2 + 1e-12;
Bi = interpolateMagneticFlux(R, X(in), Y(in), Z(in));
dsvP = [X(in) Y(in) Z(in)];  dsvB = [Bi.Bx(:) Bi.By(:) Bi.Bz(:)];
Bc = interpolateMagneticFlux(R, 0, 0, 0);  B0vec = [Bc.Bx Bc.By Bc.Bz];
lineT = linspace(-25e-3, 25e-3, 201)';
lineB = cell(1, 3);
lineP = {[lineT 0*lineT 0*lineT], [0*lineT lineT 0*lineT], [0*lineT 0*lineT lineT]};
for i = 1:3
    b = interpolateMagneticFlux(R, lineP{i}(:,1), lineP{i}(:,2), lineP{i}(:,3));
    lineB{i} = [b.Bx(:) b.By(:) b.Bz(:)];
end
planeS = linspace(-0.095, 0.095, 81);  [S1, S2] = ndgrid(planeS, planeS);
b = interpolateMagneticFlux(R, S1(:), S2(:), 0*S1(:));  planeXY = [b.Bx(:) b.By(:) b.Bz(:)];
b = interpolateMagneticFlux(R, S1(:), 0*S1(:), S2(:));  planeXZ = [b.Bx(:) b.By(:) b.Bz(:)];
save(fullfile(here, 'fem_result.mat'), 'dsvP', 'dsvB', 'B0vec', 'lineT', 'lineB', 'planeS', ...
     'planeXY', 'planeXZ', 'hDSV', 'hMag', 'hAir', 'hGrad', 'Br', 'muRec', 'dsvDiam', 'dsvLength');
fprintf('Sampling + save: %.1f s\n', toc(tSolve));
fem_postprocess;
