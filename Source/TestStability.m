alpha = 1.0;
beta = 1.0;
gamma = 1.0;
eps = 0.005;
amplitudeFn = @(i, j) alpha ./ (beta + i.^2 + j.^2) .^ (gamma/2);
grf = GRF2DON(8, amplitudeFn, [0, 0, 5, 5]);

vCoef1 = grf.generateCoef();
wCoef1 = grf.generateCoef();

voref1 = @(x, y) grf.eval(x, y, vCoef1);
woref1 = @(x, y) grf.eval(x, y, wCoef1);

[X, Y] = meshgrid(0:0.01:1.0, 0:0.01:1.0);

vCoef1_perturbed = vCoef1 + eps * randn(size(vCoef1));
voref_perturbed = @(x, y) grf.eval(x, y, vCoef1_perturbed);
figure;
surf(X, Y, voref1(X, Y), 'FaceColor', [1, 0, 0]);
hold on;
surf(X, Y, voref_perturbed(X, Y), 'FaceColor', [0, 0, 1], 'FaceAlpha', 0.2);
hold off;

fprintf("vo norm diff = %f\n", sqrt(mean(mean((voref1(X, Y) - voref_perturbed(X, Y)).^2))));

[u1, v1] = runSim(voref1, woref1);
[u2, v2] = runSim(voref_perturbed, woref1);

fprintf("final norm diff = %f\n", sqrt(mean((v1 - v2).^2)));
figure;
surf(reshape(v1, 200, 200), 'FaceColor', [1, 0, 0]);
hold on;
surf(reshape(v2, 200, 200), 'FaceColor', [0, 0, 1], 'FaceAlpha', 0.2);

function [u, v] = runSim(voref, woref)
    sim = Simulation();
    
    sim.flag_DynamicPlotting = false;
    sim.flag_video = false;
    sim.flag_ShowProgress = false;
    
    sim.bvfunc = @(x, y, t) 0.0*x;
    sim.bwfunc = @(x, y, t) 0.0*y;
    
    sim.Ti = 0.0;
    sim.Tf = 1.8;
    sim.dt = 0.01;
    sim.LoadGrid("GridFile1.mat");
    Xo = 0.0;
    Xn = 5.0;
    Yo = 0.0;
    Yn = 5.0;
    sim.ComputePDConstants();

    % x-component of initial displacement
    vofunc = @(x, y) voref((x - Xo) / (Xn - Xo), (y - Yo) / (Yn - Yo));
    % y-component of initial displacement
    wofunc = @(x, y) woref((x - Xo) / (Xn - Xo), (y - Yo) / (Yn - Yo));

    Vvofunc = @(x,y) 0.*x + 0.*y; % x-component of initial velocity
    Vwofunc = @(x,y) 0.*x + 0.*y; % y-component of initial velocity

    sim.ImposeInitialConditions(vofunc, wofunc, Vvofunc, Vwofunc);
    u = [sim.v, sim.w];  % (numNodes, 2)
    
    sim.Solve();
    v = sim.v;
end