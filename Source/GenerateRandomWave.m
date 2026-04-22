sim = Simulation();

sim.flag_DynamicPlotting = 0;
sim.flag_video = 0;
sim.flag_ShowProgress = 0;

sim.bvfunc = @(x, y) 0.0*x;
sim.bwfunc = @(x, y) 0.0*y;

sim.Ti = 0.0;
sim.Tf = 1.8;
sim.dt = 0.01;
sim.LoadGrid("GridFile1.mat");
sim.ComputePDConstants();

numChunks = 32;
datasetSize = 3000;

runSim = RandomWaveSimulator();
runSim.Setup(sim, numChunks);

xy = [sim.xx, sim.yy];
GenerateDataset("RandomWavePropagation1.ol.h5", datasetSize, numChunks, sim, xy, xy, 2, 2, @GenerateOne);

function GenerateOne(simulation, outputFile, ~, sampleIndex, curStream)
    % x-component of initial displacement
    voref = GRF2D(self.num_modes, @self.Amplitude, curStream);  
    vofunc = @(x, y) voref((x - Xo) / (Xn - Xo), (y - Yo) / (Yn - Yo));
    % y-component of initial displacement
    woref = GRF2D(self.num_modes, @self.Amplitude, curStream);
    wofunc = @(x, y) woref((x - Xo) / (Xn - Xo), (y - Yo) / (Yn - Yo));

    Vvofunc = @(x,y) 0.*x + 0.*y; % x-component of initial velocity
    Vwofunc = @(x,y) 0.*x + 0.*y; % y-component of initial velocity

    simulation.ImposeInitialConditions(vofunc, wofunc, Vvofunc, Vwofunc);
    u = [simulation.u, simulation.v];  % (numNodes, 2)
    
    simulation.Solve();
    v = [simulation.u, simulation.v];  % (numNodes, 2)
    
    WriteOLPDSample(outputFile, sampleIndex, u, v);
end