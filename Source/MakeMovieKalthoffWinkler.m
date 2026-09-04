sim = Simulation();

sim.flag_DynamicPlotting = 0;
sim.flag_video = 0;
sim.flag_ShowProgress = 1;
sim.flag_BB = true;

sim.bvfunc = @(x, y, t) 0.0*x;
sim.bwfunc = @(x, y, t) 0.0*y;

sim.Ti = 0.0;
sim.Tf = 70e-6;
sim.dt = 1e-7;
saveInterval = 25;
numFrames = 1 + floor(ceil((sim.Tf - sim.Ti) / sim.dt) / saveInterval);
fprintf("Saving %d frames\n", numFrames);
Xo = -0.1;
Xn = 0.1;
Yo = -0.05;
Yn = 0.05;
sim.del = 0.002;

sim.E = 191e9;
sim.Go = 42408;
sim.rho = 8000;
velocity = 24;
Nx = 300;
Ny = 150;
sim.ComputePDConstants();
sim.LoadGrid("KWGrid.mat");
sim.mask_nofail = zeros(size(sim.xx));
dx = (Xn - Xo) / Nx;
dy = (Yn - Yo) / Ny;
sim.CreatePreNotches([Xo/4, Yo, Xo/4, 0; Xn/4, Yo, Xn/4, 0], true);

xy = [sim.xx, sim.yy];
mask_bc = (Xo/4  - 0.25*dx < sim.xx) & (sim.xx < Xn/4 + 0.25*dx) & (Yo - dy < sim.yy) & (sim.yy < Yo + dy);

vofunc = @(x,y) 0.*x + 0.*y; % x-component of initial displacement
wofunc = @(x,y) 0.*x + 0.*y; % y-component of initial displacement
Vvofunc = @(x,y) 0.*x + 0.*y; % x-component of initial velocity
Vwofunc = @(x,y) mask_bc * velocity; % y-component of initial velocity
sim.ImposeInitialConditions(vofunc, wofunc, Vvofunc, Vwofunc);

CreateOLPDDataset("movieKW.ol.h5", xy, xy, 2, 2, numFrames);
sim.Solve(MakeCallback(mask_bc, velocity, saveInterval));

function bc = MakeCallback(mask_bc, velocity, saveInterval)
    function Callback(sim, n)
        % Enforce velocity boundary condition
        sim.Vv(mask_bc) = 0.0;
        sim.Vw(mask_bc) = velocity;
        
        if mod(n, saveInterval) == 0
            u = [sim.v, sim.w];
            v = sim.ComputeDamage();
            WriteOLPDSample("movieKW.ol.h5", 1 + round(n / saveInterval), u, v);
        end
    end
    
    bc = @Callback;
end
