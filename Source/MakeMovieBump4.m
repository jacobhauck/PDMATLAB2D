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

numFrames = 32;

xy = [sim.xx, sim.yy];

generator = MakeGenerator(numFrames);
GenerateDataset('movieBump4.ol.h5', numFrames, 1, 1234, sim, xy, xy, 2, 2, generator);


function generator = MakeGenerator(numFrames)
    function GenerateOne(simulation, outputFile, ~, sampleIndex, ~)
        % Initial displacement functions
        
        % Tolerance
        tol = 1E-15;
        
        % Parameters
        xm    = 2.5;     % x-coordinate of pulse center 
        ym    = 2.5;     % y-coordinate of pulse center 
        A     = 0.025;   % amplitude of radial Gaussian distribution
        sigma = 0.1;    % standard deviation of radial Gaussian distribution
        mu    = 6*sigma; % radial distance from pulse center of radial Gaussian distribution mean 
        
        % Functions
        r      = @(x,y) sqrt((x-xm).^2+(y-ym).^2);                    % distance from pulse center
        uo     = @(x,y) ( r(x,y) > mu-6*sigma & r(x,y) < mu+6*sigma ).* A .* exp( (-(r(x,y) - mu).^2) / (2*sigma)^2 ); % magnitude of initial displacement
        vofunc = @(x,y) uo(x,y).*(x-xm)./(r(x,y) + tol);              % x-component of initial displacement
        wofunc = @(x,y) uo(x,y).*(y-ym)./(r(x,y) + tol);              % y-component of initial displacement

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
