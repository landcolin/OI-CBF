function p = config(model, profile)
%CONFIG Parameters of the original examples and explicit paper differences.
% source: reproduce Python mathematics; paper: use literal printed formulas.
% Neither profile repairs the signed-radius issue in paper (33) silently.
if nargin < 2, profile = 'source'; end
model = validatestring(model, {'double_integrator','aircraft'});
profile = validatestring(profile, {'source','paper'});
p.model = model; p.profile = profile;
p.qpSolver = 'projection'; % Optional 'quadprog' uses Optimization Toolbox.
p.failurePolicy = 'error'; % 'backup' or explicitly 'legacy' for diagnostics
p.radiusConvention = 'printed'; % 'signed-corrected' is an explicit amendment
p.includeInitialConstraint = strcmp(profile,'paper');
p.tighteningEnabled = strcmp(profile,'source');
p.backupSaveN = 5;
if strcmp(model,'double_integrator')
    p.n = 2; p.m = 1; p.dt = .02; p.steps = floor(7/p.dt)+1;
    p.x0 = [-2;.5]; p.rtol = 1e-9; p.atol = 1e-9;
    p.T = 1; p.backupDt = .02; p.lb = -1; p.ub = 1;
    p.Lh = 1; p.supFcl = .2; p.beta = 10;
    p.alpha = @(h) 15*h+h.^3; p.alphaB = @(h) 10*h;
    p.alphaOIBackup = p.alphaB;
else
    p.n = 8; p.m = 2; p.dt = .15; p.steps = floor(110/p.dt)+1;
    p.x0 = [0;0;0;16000;0;10000;0;1]; p.rtol = 1e-6; p.atol = 1e-6;
    p.T = 60; p.backupDt = .5; p.lb = [-pi/2;-1]; p.ub = [pi/2;4];
    p.g = 9.81; p.VT = 200; p.tp = 1; p.tz = 1;
    p.kphi = 3; p.kP = 2; p.kN = .4; p.ktheta = .65;
    p.kpsi = .45; p.kH = .0001; p.Hstar = 10000;
    p.pPrimary = [27000;0]; p.Hprimary = 10000;
    p.phimax = pi/4; p.phistar = -pi/4; p.Nzstar = 1/cos(p.phistar);
    p.rho = p.VT^2/(p.g*tan(p.phistar));
    p.ng = [-1;0]; p.pg = [25000;0];
    p.cphi = pi/12; p.ctheta = pi/12; p.cH = 400;
    p.cP = pi/12; p.cNz = .1; p.c6 = 5;
    p.kappa = 100; p.beta = 1000; p.Lh = 1; p.supFcl = 150;
    p.alpha = @(h) .072*h; p.alphaB = @(h) .5*h;
    p.alphaOIBackup = p.alphaB;
    if strcmp(profile,'source'), p.alphaOIBackup = p.alpha; end
    p.altitudeNormalized = strcmp(profile,'source');
    p.absoluteHeading = strcmp(profile,'source');
end
p.backupGrid = 0:p.backupDt:p.T;
end
