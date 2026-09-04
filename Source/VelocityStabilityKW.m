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
Xo = -0.1;
Xn = 0.1;
Yo = -0.05;
Yn = 0.05;
sim.del = 0.002;

sim.E = 191e9;
sim.Go = 42408;
sim.rho = 8000;
velocityRange = [16, 32];
velocitySteps = 256;
Nx = 300;
Ny = 150;
sim.ComputePDConstants();
sim.LoadGrid("KWGrid.mat");
xy = [sim.xx, sim.yy];
sim.mask_nofail = zeros(size(sim.xx));
dx = (Xn - Xo) / Nx;
dy = (Yn - Yo) / Ny;
mask_bc = (Xo/4  - 0.25*dx < sim.xx) & (sim.xx < Xn/4 + 0.25*dx) & (Yo - dy < sim.yy) & (sim.yy < Yo + dy);

CreateOLPDDataset("velocityKW.ol.h5", xy, xy, 2, 2, numFrames);
velocities = linspace(velocityRange(1), velocityRange(2), velocitySteps);
generator = MakeGenerator(mask_bc, velocities, Xo, Xn, Yo);
GenerateDataset("velocityKW.ol.h5", velocitySteps, 1, 1234, sim, zeros(1, 1), xy, 1, 1, generator);


function generator = MakeGenerator(mask_bc, velocities, Xo, Xn, Yo)
    function GenerateOne(simulation, outputFile, ~, sampleIndex, ~)
        % Reset bonds
        simulation.LoadGrid("KWGrid.mat");  
        simulation.CreatePreNotches([Xo/4, Yo, Xo/4, 0; Xn/4, Yo, Xn/4, 0], true);
        
        % Initial conditions
        velocity = velocities(sampleIndex);
        vofunc = @(x,y) 0.*x + 0.*y; % x-component of initial displacement
        wofunc = @(x,y) 0.*x + 0.*y; % y-component of initial displacement
        Vvofunc = @(x,y) 0.*x + 0.*y; % x-component of initial velocity
        Vwofunc = @(x,y) mask_bc * velocity; % y-component of initial velocity
        simulation.ImposeInitialConditions(vofunc, wofunc, Vvofunc, Vwofunc);

        u = velocity;
        simulation.Solve(MakeCallback(mask_bc, velocity));
        v = simulation.ComputeDamage();
        
        WriteOLPDSample(outputFile, sampleIndex, u, v);
    end
    
    generator = @GenerateOne;
end

function bc = MakeCallback(mask_bc, velocity)
    function Callback(sim, ~)
        % Enforce velocity boundary condition
        sim.Vv(mask_bc) = 0.0;
        sim.Vw(mask_bc) = velocity;
    end
    
    bc = @Callback;
end
