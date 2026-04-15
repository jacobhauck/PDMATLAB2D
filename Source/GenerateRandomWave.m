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

GenerateDataset("RandomWavePropagation1.ol.h5", datasetSize, numChunks, sim, runSim);
