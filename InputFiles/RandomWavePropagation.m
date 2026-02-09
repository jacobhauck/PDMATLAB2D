
% ------------------------------------------------------------------------
%                  Domain geometry and discretization
% ------------------------------------------------------------------------

% Domain boundaries
Xo = 0; % Left  boundary of the domain
Xn = 5; % Right boundary of the domain
Yo = 0; % Lower boundary of the domain
Yn = 5; % Upper boundary of the domain

% Number of nodes
Nx = 200; % Number of nodes in the x-direction 
Ny = 200; % Number of nodes in the y-direction

% Grid perturbation coefficient
PG = 0;

GridFile = 'GridFile1.mat';

% flag_SaveGrid = 1;
% GridFileOut = 'GridFile1.mat';

% ------------------------------------------------------------------------
%                         Time discretization
% ------------------------------------------------------------------------

% Initial time
Ti = 0;

% Final time
Tf = 1.8;

% Time step
dt = 0.01;

% Time-integration scheme
TimeScheme = 'VVerlet'; % Velocity Verlet

% ------------------------------------------------------------------------
%                             PD model 
% ------------------------------------------------------------------------

% Constitutive model
model = 'GPMB'; % Generalized Prototype Microelastic Brittle (GPMB) model

% Plane elasticity model
PlanarModel = 'PlaneStress';

% Horizon
del = 0.1; 

% Influence function order indicator
omega = 0; % Constant influence function

% Flag for bond breaking
flag_BB = 0;

% ------------------------------------------------------------------------
%                      Classical material properties
% ------------------------------------------------------------------------

% Mass density
rho = 1;

% Young's modulus
E = 1;

% Fracture energy
Go = 1;

% ------------------------------------------------------------------------
%                         Meshfree discretization 
% ------------------------------------------------------------------------

% Algorithm for computation of neighbor areas
AlgName = 'FA'; % FA algorithm

% ------------------------------------------------------------------------
%                           Problem settings
% ------------------------------------------------------------------------

% ------------------
% Body force density
% ------------------

bvfunc = @(x,y,t) (0.*x + 0.*y)*t; % x-component of body force density
bwfunc = @(x,y,t) (0.*x + 0.*y)*t; % y-component of body force density

% ------------------
% Initial conditions
% ------------------ 

% Initial displacement functions

% Parameters
num_modes = 16;
alpha = 0.01;
beta = 1.0;
gamma = 1.9;
amplitude_fn = @(i, j) alpha * (beta + i.^2 + j.^2) .^ (-gamma / 2);

% Initial displacement functions
% x-component of initial displacement
voref = GRF2D(num_modes, amplitude_fn);  
vofunc = @(x, y) voref((x - Xo) / (Xn - Xo), (y - Yo) / (Yn - Yo));
% y-component of initial displacement
woref = GRF2D(num_modes, amplitude_fn);
wofunc = @(x, y) woref((x - Xo) / (Xn - Xo), (y - Yo) / (Yn - Yo));

% Initial velocity functions (zero)
Vvofunc = @(x,y) 0.*x + 0.*y; % x-component of initial velocity
Vwofunc = @(x,y) 0.*x + 0.*y; % y-component of initial velocity


% ------------------------------------------------------------------------
%                           Postprocessing
% ------------------------------------------------------------------------

% Flag for plotting during time integration
flag_DynamicPlotting = 0;

% Frequency of plotting during time integration
DynamicPlotFrequency = 10; % Plot every 10 time steps (beginning from the first one)

% Frequency of time-integration step display
TimeStepDisplayFrequency = 50000;  % don't show time-step updates

% Flag for plotting at final time
flag_FinalPlots = 0;

% Plot settings
%                     Field Name             Field variable         Colorbar title   Point size  Colormap limits   Colormap    Axes limits    Configuration
PlotSettings = {'DisplacementMagnitude' , 'sqrt(v.^2 + w.^2)/alpha'  , '$\|{\bf u}\|/\alpha$' ,    8    ,     [0.0 5.0]   ,   'parula'  , [Xo Xn Yo Yn] , 'Reference'};

% Flag to create video(s): works only if flag_DynamicPlotting = 1
flag_video = 0;

% Video frame rate
video_frate = 20;

% ------------------------------------------------------------------------