% Generates a dataset whose input u is the initial displacement field
% and output v is the displacement field at a later time (sim.Tf).
% The domain is a 5x5 square, discretized using a 200x200 grid.
% The input displacement is sampled from a Gaussian random field, and bond
% breaking is disabled.
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

num_modes = 16;
alpha = 0.01;
beta = 1.0;
gamma = 3.1;

numChunks = 64;
datasetSize = 10000;

datasetFile = "train8.ol.h5";
rngSeed = 2026;

xy = [sim.xx, sim.yy];
generator = MakeGenerator(num_modes, alpha, beta, gamma, Xo, Xn, Yo, Yn);
GenerateDataset(datasetFile, datasetSize, numChunks, rngSeed, sim, xy, xy, 2, 2, generator);


function generator = MakeGenerator(num_modes, alpha, beta, gamma, Xo, Xn, Yo, Yn)
    function GenerateOne(simulation, outputFile, ~, sampleIndex, curStream)
        amplitude = @(gx, gy) alpha * (beta + gx.^2 + gy.^2) .^ (-gamma/2);
        % x-component of initial displacement
        voref = GRF2D(num_modes, amplitude, curStream);
        vofunc = @(x, y) voref((x - Xo) / (Xn - Xo), (y - Yo) / (Yn - Yo));
        % y-component of initial displacement
        woref = GRF2D(num_modes, amplitude, curStream);
        wofunc = @(x, y) woref((x - Xo) / (Xn - Xo), (y - Yo) / (Yn - Yo));
    
        Vvofunc = @(x,y) 0.*x + 0.*y; % x-component of initial velocity
        Vwofunc = @(x,y) 0.*x + 0.*y; % y-component of initial velocity
    
        simulation.ImposeInitialConditions(vofunc, wofunc, Vvofunc, Vwofunc);
        u = [simulation.v, simulation.w];  % (numNodes, 2)
        
        simulation.Solve();
        v = [simulation.v, simulation.w];  % (numNodes, 2)
        
        WriteOLPDSample(outputFile, sampleIndex, u, v);
    end
    
    generator = @GenerateOne;
end
