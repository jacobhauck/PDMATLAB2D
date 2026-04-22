sim = Simulation();

sim.flag_DynamicPlotting = 0;
sim.flag_video = 0;
sim.flag_ShowProgress = 0;

zero = @(x, y) 0.0*x;
sim.bvfunc = zero;

sim.Ti = 0.0;
sim.Tf = 3.2e-5; % [s] :    43 microsec
sim.dt = 6.7E-8; % [s] : 67E-3 microsec
sim.del = 0.001;
sim.flag_BB = true;
sim.rho = 2440; % [kg/m^3]
sim.E = 72e+9;  % [Pa] : 72 GPa 
sim.Go = 3.8;   % [J/m^2]

sim.LoadGrid("GridFile1.mat");
PreNotchCoordinates = [-0.05  0.0  0.0  0.0];
sim.CreatePreNotches(PreNotchCoordinates, true);

Xo = min(min(sim.xx));
Xn = max(max(sim.xx));
Nbx = 129;
bx = linspace(Xo, Xn, Nbx)';  % (Nbx, 1)

Yn = max(max(abs(sim.yy)));
sim.mask_nofail = (abs(sim.yy) > Yn - sim.del); 

sim.ComputePDConstants();

numChunks = 32;
datasetSize = 3000;

alpha = 0.1;
beta = 1.0;
gamma = 1.0;

runSim = RandomWaveSimulator();
runSim.Setup(sim, numChunks);

GenerateDataset("RandomWavePropagation1.ol.h5", datasetSize, numChunks, sim, x, sim.yy, 2, 1, @GenerateOne);

function GenerateOne(simulation, outputFile, ~, sampleIndex, curStream)
    modtop_n = GRF1D(3, @(i) alpha ./ ((beta + i.^2).^(gamma/2)), curStream);
    modbot_n = GRF1D(3, @(i) alpha ./ ((beta + i.^2).^(gamma/2)), curStream);
    modtop = @(x) 6*alpha + modtop_n((x - Xo) / (Xn - Xo));
    modbot = @(x) 6*alpha + modbot_n((x - Xo) / (Xn - Xo));
    
    % y-component of body force density
    simulation.bwfunc = @(x,y,t) ( (y > Yn - dy) .* modtop(x) + (y < Yo + dy) .* modbot(x) ).*sigma.*sign(y)/dy;
    
    % Reset simulation state (need to call LoadGrid again to reset broken
    % bonds)
    simulation.ImposeInitialConditions(zero, zero, zero, zero);
    simulation.LoadGrid("GridFile1.mat");
    simulation.CreatePreNotches(PreNotchCoordinates, true);

    % Get boundary modulators
    u = [modtop_n(bx), modbot_n(bx)];  % (Nbx, 2)
    
    % Solve
    simulation.Solve();

    % Get final damage field
    v = simulation.ComputeDamage();  % (numNodes, 1)
    
    WriteOLPDSample(outputFile, sampleIndex, u, v);
end