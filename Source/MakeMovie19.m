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

numModes = 16;
beta = 1.0;
gamma = 2.5;
scale = 0.0075;
area = (Xn - Xo) * (Yn - Yo);
grf = GRF2DON(numModes, @(i, j) 1.0 ./ ((beta + i.^2 + j.^2) .^ (gamma/2)), [Xo, Yo, Xn, Yn]);
totalVar = sum(grf.amplitudeFn(grf.gx(), grf.gy()) .^ 2);
alpha = sqrt(area * scale^2 / totalVar);
fprintf("Computed alpha = %f\n", alpha);

numFrames = 32;

xy = [sim.xx, sim.yy];

generator = MakeGenerator(numModes, alpha, beta, gamma, Xo, Xn, Yo, Yn, numFrames);
GenerateDataset('movie19.ol.h5', numFrames, 1, 1234, sim, xy, xy, 2, 2, generator);


function generator = MakeGenerator(numModes, alpha, beta, gamma, Xo, Xn, Yo, Yn, numFrames)
    vofunc = [];
    wofunc = [];
    
    function GenerateOne(simulation, outputFile, ~, sampleIndex, curRngStream)
        if sampleIndex == 1
            grf = GRF2DON(numModes, @(i, j) alpha ./ ((beta + i.^2 + j.^2) .^ (gamma/2)), [Xo, Yo, Xn, Yn], curRngStream);
            % x-component of initial displacement
            vofunc = grf.generate();
            % y-component of initial displacement
            wofunc = grf.generate();
        end

        Vvofunc = @(x,y) 0.*x + 0.*y; % x-component of initial velocity
        Vwofunc = @(x,y) 0.*x + 0.*y; % y-component of initial velocity
        
        tRel = (sampleIndex - 1) / (numFrames - 1);
        simulation.Tf = simulation.Ti + (simulation.Tf - simulation.Ti) * tRel;
        simulation.ImposeInitialConditions(vofunc, wofunc, Vvofunc, Vwofunc);
        u = [simulation.v, simulation.w];  % (numNodes, 2)
        
        simulation.Solve();
        v = [simulation.v, simulation.w];  % (numNodes, 2)
        
        WriteOLPDSample(outputFile, sampleIndex, u, v);
    end
    
    generator = @GenerateOne;
end
