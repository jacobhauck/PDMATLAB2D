sim = Simulation();

sim.flag_DynamicPlotting = false;
sim.flag_video = false;
sim.flag_ShowProgress = true;

zero_t = @(x, y, t) 0.0*x;
sim.bvfunc = zero_t;

sim.Ti = 0.0;
sim.Tf = 1.0;
sim.AlgName = 'IPA-AC';
sim.flag_BB = true;

n = 100;
so = 0.05;
sim.InitDimensionless(n, so);

Xo = -0.5;
Xn =  0.5;
Yo = -0.2;
Yn =  0.2;

Nx = 300;
Ny = 120;

if exist("ConvergenceGridDimensionless.mat", "file")
    sim.LoadGrid("ConvergenceGridDimensionless.mat");
else
    sim.GenerateGrid(Xo, Xn, Yo, Yn, Nx, Ny, false);
    sim.SaveGrid("ConvergenceGridDimensionless.mat");
end

PreNotchCoordinates = [-0.5  0.0  0.0  0.0];
sim.CreatePreNotches(PreNotchCoordinates, true);

Nbx = 33;
bx = linspace(Xo, Xn, Nbx)';  % (Nbx, 1)

sim.mask_nofail = (abs(sim.yy) > Yn - sim.del); 

sim.ComputePDConstants();

numChunks = 1;
seed = 1234;

sigma = 0.00025;
sigmaTop = 0.3;
sigmaBot = 0.6;
dt = [1/200, 1/400];
datasetSize = length(dt);

datasetName = "convergenceDimensionless.ol.h5";

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
        simulation.LoadGrid("ConvergenceGridDimensionless.mat");
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
