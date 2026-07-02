% Simulate wave propagation
sim = Simulation();

% Turn off progress reporting and video output
sim.flag_DynamicPlotting = 0;
sim.flag_video = 0;
sim.flag_ShowProgress = 0;

% Set zero body force
sim.bvfunc = @(x, y, t) 0.0*x;
sim.bwfunc = @(x, y, t) 0.0*y;

% Material properties
sim.E_1 = 1.0;
sim.E_2 = 2.0;
sim.eta_12_11 = 0.0;
sim.eta_12_22 = 0.0;
sim.nu_12 = 0.1;
sim.rho = 1.0;
sim.model = 'Anisotropic';
sim.PlanarModel = 'Pure';

% Domain
sim.Ti = 0.0;
sim.Tf = 1.8;
sim.dt = 0.01;
Xo = 0.0;
Xn = 5.0;
Yo = 0.0;
Yn = 5.0;
sim.ComputePDConstants();
if ~exist("GridFile1-rectangular.mat", "file")
    sim.GenerateGrid(Xo, Xn, Yo, Yn, 200, 200, false);
    sim.SaveGrid("GridFile1-rectangular.mat");
else
    sim.LoadGrid("GridFile1-rectangular.mat");
end

% Dataset
numChunks = 1;
datasetSize = 1;
randomSeed = 1234;
datasetName = "waveRectangular.ol.h5";

xy = [sim.xx, sim.yy];
generator = MakeGenerator();
GenerateDataset(datasetName, datasetSize, numChunks, randomSeed, sim, xy, xy, 2, 2, generator);


function generator = MakeGenerator()
    function GenerateOne(simulation, outputFile, ~, sampleIndex, ~)
        % Parameters
        xm    = 2.5;     % x-coordinate of pulse center 
        ym    = 2.5;     % y-coordinate of pulse center 
        A     = 0.025;   % amplitude of radial Gaussian distribution
        sigma = 0.1;     % standard deviation of radial Gaussian distribution
        mu    = 6*sigma; % radial distance from pulse center of radial Gaussian distribution mean 
        tol   = 1e-10;
        
        % Functions
        r      = @(x,y) sqrt((x-xm).^2+(y-ym).^2);                    % distance from pulse center
        uo     = @(x,y) ( r(x,y) > mu-6*sigma & r(x,y) < mu+6*sigma ).* A .* exp( (-(r(x,y) - mu).^2) / (2*sigma)^2 ); % magnitude of initial displacement
        vofunc = @(x,y) uo(x,y).*(x-xm)./(r(x,y) + tol);              % x-component of initial displacement
        wofunc = @(x,y) uo(x,y).*(y-ym)./(r(x,y) + tol);              % y-component of initial displacement

        
        % Zero initial velocity
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
