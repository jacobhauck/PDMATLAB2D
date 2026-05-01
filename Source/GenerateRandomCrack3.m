% Generates a dataset whose input u comprises a pair of modulation
% functions used to modulate the traction on the top and bottom of a bar,
% and whose output v is the final damage field. The modulation functions
% are constant functions whose value is sampled uniformly from a range.

sim = Simulation();

sim.flag_DynamicPlotting = false;
sim.flag_video = false;
sim.flag_ShowProgress = false;

zero_t = @(x, y, t) 0.0*x;
sim.bvfunc = zero_t;

sim.Ti = 0.0;
sim.Tf = 4.3e-5; % [s] :    43 microsec
sim.dt = 6.7E-8; % [s] : 67E-3 microsec
sim.del = 0.001;
sim.flag_BB = true;
sim.rho = 2440; % [kg/m^3]
sim.E = 72e+9;  % [Pa] : 72 GPa 
sim.Go = 3.8;   % [J/m^2]

Xo = -0.05; % [m] : Left  boundary of the domain 
Xn =  0.05; % [m] : Right boundary of the domain
Yo = -0.02; % [m] : Lower boundary of the domain
Yn =  0.02; % [m] : Upper boundary of the domain

Nx = 300;
Ny = 120;
dy = (Yn - Yo) / Ny;

sim.LoadGrid("GridFile2.mat");
PreNotchCoordinates = [-0.05  0.0  0.0  0.0];
sim.CreatePreNotches(PreNotchCoordinates, true);

Nbx = 129;
bx = linspace(Xo, Xn, Nbx)';  % (Nbx, 1)

sim.mask_nofail = (abs(sim.yy) > Yn - sim.del); 

sim.ComputePDConstants();

numChunks = 64;
datasetSize = 500;
seed = 1234;

sigma = 2E6;
sigmaTopRange = [0.0, 0.38];
sigmaBotRange = [0.75, 0.75];
numModes = 3;
alpha = 1.0;
beta = 1.0;
gamma = 2.0;

datasetName = "test3.ol.h5";

xy = [sim.xx, sim.yy];
generator = MakeGenerator(numModes, alpha, beta, gamma, sigmaTopRange, sigmaBotRange, Yo, Yn, dy, sigma, bx, PreNotchCoordinates);
GenerateDataset(datasetName, datasetSize, numChunks, seed, sim, bx, xy, 2, 1, generator);

function normalized = linearNormalize(vals)
    m = min(vals);
    M = max(vals);
    normalized = (vals - m) / (M - m);
end

function generator = MakeGenerator(numModes, alpha, beta, gamma, sigmaTopRange, sigmaBotRange, Yo, Yn, dy, sigma, bx, PreNotchCoordinates)
    function GenerateOne(simulation, outputFile, ~, sampleIndex, curStream)
        amplitudeFn = @(g) alpha ./ (beta + g.^2) .^ (gamma/2);
        modtop_n = GRF1D(numModes, amplitudeFn, curStream);
        modbot_n = GRF1D(numModes, amplitudeFn, curStream);
        modtop = @(x) linearNormalize(modtop_n(x)) * (sigmaTopRange(2) - sigmaTopRange(1)) + sigmaTopRange(1);
        modbot = @(x) linearNormalize(modbot_n(x)) * (sigmaBotRange(2) - sigmaBotRange(1)) + sigmaBotRange(1);
        
        % y-component of body force density
        simulation.bwfunc = @(x,y,t) ( (y > Yn - dy) .* modtop(x) + (y < Yo + dy) .* modbot(x) ).*sigma.*sign(y)/dy;
        
        % Reset simulation state (need to call LoadGrid again to reset broken
        % bonds)
        zero = @(x, y) 0.0 * x;
        simulation.ImposeInitialConditions(zero, zero, zero, zero);
        simulation.LoadGrid("GridFile2.mat");
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
