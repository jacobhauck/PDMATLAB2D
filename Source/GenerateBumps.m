% Generates a dataset whose input u is the initial displacement field
% and output v is the displacement field at a later time (sim.Tf).
% The domain is a 5x5 square, discretized using a 200x200 grid.
% The input displacement is sampled from a Gaussian random field, and bond
% breaking is disabled.
sim = Simulation();

sim.flag_DynamicPlotting = 0;
sim.flag_video = 0;
sim.flag_ShowProgress = 0;

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

numChunks = 1;
datasetSize = 5;
datasetName = "smooth_bumps.ol.h5";
randomSeed = 1234;

bumpX = 2.5;
bumpY = 2.5;
bumpWidth = 1.5;
bumpHeight = 0.015;

xy = [sim.xx, sim.yy];
generator = MakeGenerator(bumpX, bumpY, bumpWidth, bumpHeight);
GenerateDataset(datasetName, datasetSize, numChunks, randomSeed, sim, xy, xy, 2, 2, generator);


function generator = MakeGenerator(bumpX, bumpY, width, height)
    function GenerateOne(simulation, outputFile, ~, sampleIndex, ~)
        bump = SmoothBump1d(sampleIndex - 1, width, height);
        % x-component of initial displacement
        vofunc = @(x,y) bump.evaluate(sqrt((x - bumpX).^2 + (y - bumpY).^2));
        % y-component of initial displacement
        wofunc = @(x,y) 0.*x;
    
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
