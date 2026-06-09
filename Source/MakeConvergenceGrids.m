sim = Simulation();

sim.Ti = 0.0;
sim.Tf = 4.3e-5; % [s] :    43 microsec
sim.del = 0.001;
sim.flag_BB = true;
sim.rho = 2440; % [kg/m^3]
sim.E = 72e+9;  % [Pa] : 72 GPa 
sim.Go = 3.8;   % [J/m^2]

Xo = -0.05; % [m] : Left  boundary of the domain 
Xn =  0.05; % [m] : Right boundary of the domain
Yo = -0.02; % [m] : Lower boundary of the domain
Yn =  0.02; % [m] : Upper boundary of the domain

rows = '5678';
Nx0 = 300;
dx = [(Xn - Xo) / Nx0, (Xn - Xo) / (2 * Nx0), (Xn - Xo) / (4 * Nx0), (Xn - Xo) / (8 * Nx0)];
delta = [1, 1/2, 1/4, 1/8] * 0.001;

for iRow = 1:length(rows)
    row = rows(iRow);
    Nx = round((Xn - Xo) / dx(iRow));
    Ny = round(2 * Nx / 5);
    sim.del = delta(iRow);

    fprintf('Making grid file "ConvergenceGrid%s-FA.mat" with Nx=%d, Ny=%d\n', row, Nx, Ny);
    sim.AlgName = 'FA';
    sim.GenerateGrid(Xo, Xn, Yo, Yn, Nx, Ny, false);
    sim.SaveGrid(sprintf("ConvergenceGrid%s-FA.mat", row));

    fprintf('Making grid file "ConvergenceGrid%s-IPA-AC.mat" with Nx=%d, Ny=%d\n', row, Nx, Ny);
    sim.AlgName = 'IPA-AC';
    sim.GenerateGrid(Xo, Xn, Yo, Yn, Nx, Ny, false);
    sim.SaveGrid(sprintf("ConvergenceGrid%s-IPA-AC.mat", row));
end
