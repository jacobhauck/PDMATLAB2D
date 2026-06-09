sim = Simulation();

sim.flag_DynamicPlotting = false;
sim.flag_video = false;
sim.flag_ShowProgress = false;

zero_t = @(x, y, t) 0.0*x;
sim.bvfunc = zero_t;

sim.Ti = 0.0;
sim.Tf = 4.3e-5; % [s] :    43 microsec
sim.AlgName = 'IPA-AC';
sim.del = 0.001 / 8;
sim.flag_BB = true;
sim.rho = 2440; % [kg/m^3]
sim.E = 72e+9;  % [Pa] : 72 GPa 
sim.Go = 3.8;   % [J/m^2]

Xo = -0.05; % [m] : Left  boundary of the domain 
Xn =  0.05; % [m] : Right boundary of the domain
Yo = -0.02; % [m] : Lower boundary of the domain
Yn =  0.02; % [m] : Upper boundary of the domain

Nx = 300*8;
Ny = 120*8;

sim.LoadGrid("ConvergenceGrid8-IPA-AC.mat");
PreNotchCoordinates = [-0.05  0.0  0.0  0.0];
sim.CreatePreNotches(PreNotchCoordinates, true);

Nbx = 33;
bx = linspace(Xo, Xn, Nbx)';  % (Nbx, 1)

sim.mask_nofail = (abs(sim.yy) > Yn - sim.del); 

sim.ComputePDConstants();

numChunks = 4;
datasetSize = 4;
seed = 1234;

sigma = 2E6;
sigmaTop = 0.32;
sigmaBot = 0.75;
dt = [6.70E-08, 3.35E-08, 1.68E-08, 8.38E-09];

datasetName = "convergence8-IPA-AC.ol.h5";

xy = [sim.xx, sim.yy];
generator = MakeGenerator(datasetSize, numChunks, dt, sigmaTop, sigmaBot, Yo, Yn, sigma, bx, PreNotchCoordinates);
GenerateDataset(datasetName, datasetSize, numChunks, seed, sim, bx, xy, 2, 1, generator);

function generator = MakeGenerator(datasetSize, numChunks, dt, sigmaTop, sigmaBot, Yo, Yn, sigma, bx, PreNotchCoordinates)
    function GenerateOne(simulation, outputFile, chunkIndex, sampleIndex, ~)
        minChunkSize = floor(datasetSize / numChunks);
        numExtras = mod(datasetSize, numChunks);
        globalIndex = (chunkIndex - 1) * minChunkSize + min(numExtras, chunkIndex - 1) + sampleIndex;
        simulation.dt = dt(globalIndex);

        modtop = @(x) ones(size(x)) * sigmaTop;
        modbot = @(x) ones(size(x)) * sigmaBot;
        
        % y-component of body force density
        simulation.bwfunc = @(x,y,t) ( (y > Yn - simulation.del) .* modtop(x) + (y < Yo + simulation.del) .* modbot(x) ).*sigma.*sign(y)/simulation.del;
        
        % Reset simulation state (need to call LoadGrid again to reset broken
        % bonds)
        zero = @(x, y) 0.0 * x;
        simulation.ImposeInitialConditions(zero, zero, zero, zero);
        simulation.LoadGrid("ConvergenceGrid8-IPA-AC.mat");
        simulation.CreatePreNotches(PreNotchCoordinates, true);
    
        % Get boundary modulators
        u = [modtop(bx), modbot(bx)];  % (Nbx, 2)
        
        % Solve
        simulation.Solve();
    
        % Get final damage field
        v = simulation.ComputeDamage();  % (numNodes, 1)
        
        WriteOLPDSample(outputFile, sampleIndex, u, v);
    end

    generator = @GenerateOne;
end
