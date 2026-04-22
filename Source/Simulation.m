classdef Simulation < handle & matlab.mixin.Copyable
    % Defines a single PDMATLAB2D simulation
    % 
    % This is an encapsulation of the global variable approach used
    % in the original implementation. Rather than writing an InputDeck,
    % one may instead construct a Simulation object containing all the
    % relevant input data. The Simulation object will also contain the
    % simulation state, which is modified using the original functional
    % interface.

    properties
        % ==================== Input Variables ====================
        
        % ===== Visualization =====

        % Whether to use dynamic plotting
        flag_DynamicPlotting = true;

        % Whether to make and output a video file
        flag_video = true;
        
        % Whether to show progress messages
        flag_ShowProgress = true;

        % Settings for generating plots
        % Cell array with columns [Field Name, Field variable, Colorbar
        % title, Point size, Colormap limits, Colormap, Axes limits,
        % Configuration]
        PlotSettings = {}

        % How often (in # of time steps) to display time step update
        % messages
        TimeStepDisplayFrequency = 10

        % How often (in # of time steps) to update dynamic plots
        DynamicPlotFrequency = 10

        % Name of simulation for output files
        OutputName

        % Video files being written
        vidfile

        % ===== Peridynamic Model =====
        
        % Constitutive model
        model = 'GPMB'; % Generalized Prototype Microelastic Brittle (GPMB) model

        % Plane elasticity model
        PlanarModel = 'PlaneStress';

        % Horizon
        del = 0.1; 

        % Influence function order indicator
        omega = 0; % Constant influence function
        
        % Flag for bond breaking
        flag_BB = false;

        % ===== Material properties =====
        
        % Mass density
        rho = 1;
        
        % Young's modulus
        E = 1;
        
        % Fracture energy
        Go = 1;

        % ===== Body forces =====
        
        % x body force density function
        bvfunc

        % y body force density function
        bwfunc

        % ===== Discretization =====

        % Algorithm to use when generating neighbor information
        AlgName = 'FA';

        % Time stepping scheme
        TimeScheme = 'VVerlet';

        % ==================== Simulation state ====================

        % ===== Discretization =====
        
        % Initial time
        Ti

        % Final time
        Tf

        % Time step
        dt

        % x coordinates, shape (N, 1) where N = total number of grid points
        xx
        
        % y coordinates, shape (N, 1) where N = total number of grid points
        yy
        
        % array of neighbor numbers for all nodes (2D array for most cases)
        u_NA

        % array of influence function values of neighbor bonds for all 
        % nodes (2D array for most cases)
        IF_NA

        % array of neighbor areas for all nodes (2D array for most cases)
        V_NA

        % array of reference lengths of neighbor bonds for all nodes 
        % (2D array for most cases)
        r_hat_NA

        % array of x-coordinates of quadrature points for all nodes 
        % (2D array for most cases)
        x_hat_NA

        % array of y-coordinates of quadrature points for all nodes 
        % (2D array for most cases)
        y_hat_NA

        % Whether a rectangular domain with uniform grid is being used 
        flag_RDUG = 0
        
        % Mask for no-fail condition using current discretization
        mask_nofail

        % Cached denominator for damage function
        phiD

        % ===== Peridynamic constants =====

        % Micromodulus constant
        c

        % Critical stretch constant
        so

        % ===== Kinematic state =====
        
        % current time step index
        k = 0

        % current time
        t

        % x displacement (on current time step)
        v

        % y displacement (on current time step)
        w

        % x velocity (on current time step)
        Vv

        % y velocity (on current time step)
        Vw
        
        % ===== Force and energy densities =====
        
        % x internal force density
        Fv

        % y internal force density
        Fw

        % internal energy density
        W

        % x body force density
        bv

        % y body force density
        bw

    end

    methods
        function self = GenerateGrid(self, Xo, Xn, Yo, Yn, Nx, Ny, PG)
            % Generate a grid
            if self.flag_ShowProgress
                tic
            end
            [self.xx, self.yy, ~, ~, dx, dy, VV, xx1, yy1, M] = GridGenerator(Xo, Xn, Yo, Yn, Nx, Ny, PG);
            if self.flag_ShowProgress
                fprintf('Generate grid ................................... = %f (sec) \n', toc);
            end
        
            % Tolerance
            tol = 1E-15;
        
            % Flag for Rectangular Domain Uniform Grid (RDUG)
            if abs(dx - dy) < tol 
                self.flag_RDUG = 1;
            else
                self.flag_RDUG = 0;
            end
        
            % ---------------------------------
            %     Generate neighbor list
            % ---------------------------------
            
            if self.flag_ShowProgress
                tic
            end
            
            [nl{1:6}] = NeighborList( ...
                Nx, Ny, ...
                self.xx, self.yy, xx1, yy1, M, ...
                self.del, ...
                dx, dy, ...
                VV, self.omega, ...
                self.AlgName, self.flag_RDUG ...
            );
            self.u_NA = nl{1};
            self.IF_NA = nl{2}; 
            self.V_NA = nl{3};
            self.r_hat_NA = nl{4};
            self.x_hat_NA = nl{5}; 
            self.y_hat_NA = nl{6};
            
            if self.flag_ShowProgress
                fprintf('Generate neighbor list .......................... = %f (sec)\n',toc);
            end
        end

        function SaveGrid(self, GridFileOut)
            save(GridFileOut, "-struct", "self", ...
                "xx", "yy", ...
                "u_NA", "IF_NA", "V_NA", ...
                "r_hat_NA", "x_hat_NA", "y_hat_NA", ...
                "flag_RDUG" ...
            );

            if self.flag_ShowProgress
                fprinf("Saved grid to file: %s\n", GridFileOut);
            end
        end

        function LoadGrid(self, GridFileIn)
            loadVars = {
                "xx", "yy", ...
                "u_NA", "IF_NA", "V_NA", ...
                "r_hat_NA", "x_hat_NA", "y_hat_NA", ...
                "flag_RDUG"
            };

            in = load(GridFileIn, loadVars{:});
            self.xx = in.xx;
            self.yy = in.yy;
            self.u_NA = in.u_NA;
            self.IF_NA = in.IF_NA;
            self.V_NA = in.V_NA;
            self.r_hat_NA = in.r_hat_NA;
            self.x_hat_NA = in.x_hat_NA;
            self.y_hat_NA = in.y_hat_NA;
            self.flag_RDUG = in.flag_RDUG;

            if self.flag_ShowProgress
                fprintf("Loaded grid from file: %s\n", GridFileIn);
            end
        end

        function CreatePreNotches(self, PreNotchCoordinates, flag_DamagedPrenotches)
            if self.flag_ShowProgress
                tic
            end
        
            % Find number of pre-notches
            [s1, ~] = size(PreNotchCoordinates);
        
            % Loop over pre-notches
            for n = 1:s1
                % Coordinates of one endpoint of the pre-notch
                Xc1 = PreNotchCoordinates(n, 1);
                Yc1 = PreNotchCoordinates(n, 2);
        
                % Coordinates of the other endpoint of the pre-notch
                Xc2 = PreNotchCoordinates(n, 3);
                Yc2 = PreNotchCoordinates(n, 4);
        
                % Create pre-notch
                [self.u_NA] = PreNotch(self.xx, self.yy, self.u_NA, Xc1, Yc1, Xc2, Yc2);
            end

            if nargin > 2
                % Compute damage ratio denominator before application of pre-notches
                if flag_DamagedPrenotches
                    self.phiD = sum(self.V_NA, 2);
                
                % Compute damage ratio denominator after application of pre-notches
                else
                    self.phiD = sum((self.u_NA>0) .* self.V_NA, 2);
                end
            else
                % Compute default damage ratio denominator (after application of pre-notches)
                self.phiD = sum((self.u_NA>0) .* self.V_NA, 2);
            end

            if self.flag_ShowProgress
                fprintf('Create pre-notch(es) ............................ = %f (sec)\n', toc)
            end
        end

        function ComputePDConstants(self)
            if self.flag_ShowProgress
                tic
            end
            
            [self.c, self.so] = PDBondConstants(self.omega, self.del, self.E, self.Go, self.model, self.PlanarModel);
            
            if self.flag_ShowProgress
                fprintf('Compute PD constants ............................ = %f (sec)\n', toc)
            end
        end

        function ImposeInitialConditions(self, vofunc, wofunc, Vvofunc, Vwofunc)
            if self.flag_ShowProgress
                tic
            end
            
            % Compute initial displacement for all nodes
            self.v = vofunc(self.xx, self.yy);   % x-component of initial displacement
            self.w = wofunc(self.xx, self.yy);   % y-component of initial displacement
            
            % Compute initial velocity for all nodes
            self.Vv = Vvofunc(self.xx, self.yy); % x-component of initial velocity
            self.Vw = Vwofunc(self.xx, self.yy); % y-component of initial velocity
            
            if self.flag_ShowProgress
                fprintf('Impose initial conditions ....................... = %f (sec)\n', toc)
            end
        end

        function ComputeForceEnergyDensity(self)
            if self.flag_ShowProgress
                tic
            end
            
            [self.Fv, self.Fw, self.W] = ForceEnergyDensity( ...
                self.xx, self.yy, ...
                self.v, self.w, ...
                self.c, self.u_NA, self.IF_NA, self.V_NA, ...
                self.r_hat_NA, self.x_hat_NA, self.y_hat_NA, ...
                self.model, self.flag_RDUG ...
            );
            
            if self.flag_ShowProgress
                fprintf('Compute internal force & energy densities = %f (sec)\n', toc)
            end
        end

        function ComputeBodyForceDensity(self)
            if isempty(self.t)
                error("Cannot compute body force density until t is set.");
            end

            if self.flag_ShowProgress
                tic
            end
            
            self.bv = self.bvfunc(self.xx, self.yy, self.t); % x-component of initial body force density
            self.bw = self.bwfunc(self.xx, self.yy, self.t); % y-component of initial body force density
            
            if self.flag_ShowProgress
                fprintf('Compute body force density .............. = %f (sec)\n', toc)
            end
        end

        function phi = ComputeDamage(self)
            phi = 1 - sum((self.u_NA>0) .* self.V_NA, 2) ./ self.phiD;
        end
        
        function CreateVideoFiles(self)
            % Check if video output directory exists and create it otherwise
            if ~exist(['../Videos/' self.OutputName], 'dir')
                mkdir(['../Videos/' self.OutputName])
            end

            % Find number of plots: a video is generated for each plot
            [s1,~] = size(self.PlotSettings);

            % Initialize VideoWriter array
            self.vidfile = VideoWriter.empty(s1, 0);

            % Loop over plots
            for nplot = 1:s1
                % Read plot field name
                field_name = self.PlotSettings{nplot, 1};

                % Video file name
                video_filename = ['../Videos/' self.OutputName '/' field_name '.mp4'];

                % Open video file
                self.vidfile(nplot) = VideoWriter(video_filename, 'MPEG-4');
                self.vidfile(nplot).FrameRate = video_frate;
                open(self.vidfile(nplot));
            end

        end

        function PlotFields(self, flag_Save)
            % Find number of plots
            [s1,~] = size(self.PlotSettings);

            % Loop over plots
            for nplot = 1:s1
                % Figure number
                nfig = nplot;                       
                
                % Read plot field variable
                cnodes = eval(self.PlotSettings{nplot, 2}); 

                % Read plot settings
                ctitle        = self.PlotSettings{nplot, 3};       % Colorbar title
                psize         = self.PlotSettings{nplot, 4};       % Point size
                climits       = self.PlotSettings{nplot, 5};       % Colormap limits
                cmap          = self.PlotSettings{nplot, 6};       % Colormap
                box           = self.PlotSettings{nplot, 7};       % Axes limits
                configuration = self.PlotSettings{nplot, 8};       % Configuration: 'Reference' or 'Current'

                % Plot field
                if strcmp(configuration,'Reference')
                    PlotField(nfig, self.xx, self.yy, cnodes, ctitle, psize, climits, cmap, box)
                elseif strcmp(configuration,'Current')
                    PlotField(nfig, self.xx + self.v, self.yy + self.w, cnodes, ctitle, psize, climits, cmap, box)
                else
                    error('Invalid configuration.');
                end

                if flag_Save
                    % Save figure
                    if ~exist(['../Outputs/' self.OutputName '/' cnodes_name], 'dir')
                        mkdir(['../Outputs/' self.OutputName '/' cnodes_name])
                    end
                    filename = ['../Outputs/' self.OutputName '/' cnodes_name '.eps'];
                    saveas(gcf, filename, 'epsc');
                end

                % Video update
                if self.flag_video && self.flag_DynamicPlotting
                    % Access figure window
                    figure(nfig);

                    % Format time for title
                    time = sprintf('%.2e', self.t);
                    time = regexprep(time, '(e[\+\-])0(\d)', '$1$2');

                    % Set figure title
                    title(['Time = ' time '\,s'], 'Interpreter', 'latex')

                    % Capture frame and write it to video file
                    frame = getframe(gcf);
                    writeVideo(self.vidfile(nplot), frame);
                    
                    if flag_Save
                        close(self.vidfile(nplot))
                    end
                end
            end
        end

        function Solve(self)
            self.k = 0;
            self.t = self.Ti;
            tVec = self.Ti : self.dt : self.Tf;
            Nt = length(tVec);
            
            if self.flag_ShowProgress
                tic;
            end
            
            self.ComputeForceEnergyDensity();
            self.ComputeBodyForceDensity();

            % Loop over time steps
            for n = 1:Nt-1
                [nextstate{1:10}] = TimeIntegrator( ...
                    self.TimeScheme, ...
                    self.xx, self.yy, ...
                    self.v, self.w, ...
                    self.Vv, self.Vw, ...
                    self.Fv, self.Fw, ...
                    self.bv, self.bw, ...
                    self.t, self.bvfunc, self.bwfunc, ...
                    self.dt, ...
                    self.u_NA, self.IF_NA, self.V_NA, ...
                    self.r_hat_NA, self.x_hat_NA, self.y_hat_NA, ...
                    self.rho, self.c, self.model, ...
                    self.flag_RDUG, self.so, ...
                    self.mask_nofail, self.flag_BB ...
                );

                self.v = nextstate{1};
                self.w = nextstate{2};
                self.Vv = nextstate{3};
                self.Vw = nextstate{4};
                self.Fv = nextstate{5};
                self.Fw = nextstate{6};
                self.bv = nextstate{7};
                self.bw = nextstate{8};
                self.W = nextstate{9};
                self.u_NA = nextstate{10};
            
                % Time-integration step display
                if self.flag_ShowProgress && (mod(n, self.TimeStepDisplayFrequency) == 0)
                    fprintf('Time integration (n = %4g / %g ) = %f (sec)\n', n, Nt-1, toc)
                end
            
                % Time after performing time-integration step
                self.t = tVec(n+1);
                self.k = self.k + 1;
            
                if self.flag_DynamicPlotting && (n == 1 || mod(n-1, self.DynamicPlotFrequency) == 0)
                    self.PlotFields(false);
                end
            end

            if self.flag_DynamicPlotting
                self.PlotFields(true);
            end
        end
    end
end