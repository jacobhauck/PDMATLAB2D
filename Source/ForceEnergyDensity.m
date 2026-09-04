
% ========================================================================
% Copyright (c) 2022 by Oak Ridge National Laboratory                      
% All rights reserved.                                                     
%                                                                           
% This file is part of PDMATLAB2D. PDMATLAB2D is distributed under a           
% BSD 3-clause license. For the licensing terms see the LICENSE file in    
% the top-level directory.                                                 
%                                                                          
% SPDX-License-Identifier: BSD-3-Clause                                    
% ========================================================================

% ========================================================================
% The function ForceEnergyDensity computes the internal force density and 
% the macroelastic energy density for all nodes
% ========================================================================

% Input
% -----
% xx         : x coordinates of all nodes in the grid 
% yy         : y coordinates of all nodes in the grid 
% v          : displacement of each node in the x-direction
% w          : displacement of each node in the y-direction
% c          : micromodulus constant (or array of micromodulus function
%              per-bond, for Anisotropic model)
% u_NA       : array of neighbor numbers for all nodes
% IF_NA      : array of influence function values of neighbor bonds for all nodes 
% V_NA       : array of neighbor areas for all nodes 
% r_hat_NA   : array of reference lengths of neighbor bonds for all nodes 
% x_hat_NA   : array of x-coordinates of quadrature points for all nodes
% y_hat_NA   : array of y-coordinates of quadrature points for all nodes 
% model      : constitutive model
%              'GPMB'         : generalized prototype microelastic brittle
%              'Anisotropoic' : anisotropic PMB model
%              'Classical'    : classical elasticity model with homogeneous
%                               Dirichlet boundary conditions
% flag_RDUG  : if == 1, then the grid is uniform over a rectangular domain
%              (dx = dy); the flag name RDUG stands for "Rectangular Domain
%              Uniform Grid" 
%              if == 0, then the grid is a general grid
% Nx         : Number of nodes in x direction. Only used for 'Classical'
%              model
% E          : Young's modulus, if using 'Classical' model

% Output
% ------
% Fv         : x-component of internal force density for all nodes
% Fw         : y-component of internal force density for all nodes
% W          : macroelastic energy density for all nodes

% Discussion:
% ----------
% The GPMB model is a generalization of the PMB model by incorporating an 
% influence function. 
% 
% The PMB model was presented in:
%
% S.A. Silling and E. Askari, A meshfree method based on the peridynamic 
% model of solid mechanics, Computers and Structures 83 (2005): 1526–1535.
% 
% The GPMB model was presented in:
%
% P. Seleson and M. L. Parks, On the role of the influence function in the 
% peridynamic theory, International Journal for Multiscale Computational 
% Engineering 9(6) (2011): 689–706.
%
%

function [Fv,Fw,W] = ForceEnergyDensity(xx,yy,v,w,c,u_NA,IF_NA,V_NA,r_hat_NA,x_hat_NA,y_hat_NA,model,flag_RDUG,Nx,E)

    % --------------------------------------------------------------------
    %                           GPMB model
    % --------------------------------------------------------------------

    if strcmp(model,'GPMB')
    
        % Number of nodes
        Nnodes = length(xx);
    
        % Initialize internal force density components and macroelastic energy density arrays
        Fv = zeros(Nnodes,1); % Array of x-components of internal force density for all nodes
        Fw = zeros(Nnodes,1); % Array of y-components of internal force density for all nodes
        W  = zeros(Nnodes,1); % Array of macroelastic energy density for all nodes
   
        % Find maximum number of neighbors any node can have
        zmax = length(u_NA(1,:));
    
        if flag_RDUG == 1

            % --------------------------------------------------------
            %              Computation for uniform grid
            %                 over rectangular domain
            % --------------------------------------------------------

            % Loop through all nodes ui
            for ui = 1:Nnodes

                % Get x-coordinate of node ui
                xi = xx(ui);
                % Get y-coordinate of node ui
                yi = yy(ui);

                % Get x-component of displacement of node ui
                vi = v(ui);
                % Get y-component of displacement of node ui
                wi = w(ui);

                % Loop over neighbors of node ui
                for z = 1:zmax

                    % Get neighbor cell number
                    uk = u_NA(ui,z);

                    % Only consider neighbor cells with a larger number
                    % than the source node to compute the bond pairwise force
                    % and pairwise potential once per bond
                    if uk > ui

                        % Get x-coordinate of quadrature point
                        xk_hat = x_hat_NA(ui,z);
                        % Get y-coordinate of quadrature point
                        yk_hat = y_hat_NA(ui,z);

                        % Get x-component of displacement of node uk
                        vk = v(uk);
                        % Get y-component of displacement of node uk
                        wk = w(uk);

                        % Read per-bond quantities for bond ui-uk
                        IFk = IF_NA(ui,z);         % Influence function
                        Vk  = V_NA(ui,z);          % Neighbor area
                        Rk  = r_hat_NA(ui,z);      % Reference bond length

                        % Compute deformed bond components and length
                        rxk = xk_hat-xi + vk-vi;   % x-component of current relative position
                        ryk = yk_hat-yi + wk-wi;   % y-component of current relative position
                        rk  = sqrt(rxk^2 + ryk^2); % current bond length

                        % Compute bond stretch
                        sk = (rk - Rk)/Rk;

                        % ------------------------------------------------
                        %        Compute internal force density
                        % ------------------------------------------------

                        % Compute bond pairwise force magnitude
                        fk = c*IFk*sk;

                        % Compute bond pairwise force components
                        fvk = fk*rxk/rk;
                        fwk = fk*ryk/rk;

                        % Update internal force function components of node ui
                        Fv(ui) = Fv(ui) + fvk*Vk;
                        Fw(ui) = Fw(ui) + fwk*Vk;

                        % Update internal force function components of node uk
                        Fv(uk) = Fv(uk) - fvk*Vk;
                        Fw(uk) = Fw(uk) - fwk*Vk;

                        % ------------------------------------------------
                        %      Compute macroelastic energy density
                        % ------------------------------------------------

                        % Compute bond pairwise potential
                        wk = 0.5*c*IFk*(sk^2)*Rk;

                        % Update macroelastic energy density of node ui
                        W(ui) = W(ui) + 0.5*wk*Vk;

                        % Update macroelastic energy density of node uk
                        W(uk) = W(uk) + 0.5*wk*Vk;

                    end
                end
            end

        elseif flag_RDUG == 0    

            % --------------------------------------------------------
            %              Computation for general grids
            % --------------------------------------------------------

            % Loop through all nodes ui
            for ui = 1:Nnodes

                % Get x-coordinate of node ui
                xi = xx(ui);
                % Get y-coordinate of node ui
                yi = yy(ui);

                % Get x-component of displacement of node ui
                vi = v(ui);
                % Get y-component of displacement of node ui
                wi = w(ui);

                % Loop over neighbors of node ui
                for z = 1:zmax

                    % Get neighbor cell number
                    uk = u_NA(ui,z);

                    % Only consider neighbor cells
                    if uk > 0

                        % Get x-coordinate of quadrature point
                        xk_hat = x_hat_NA(ui,z);
                        % Get y-coordinate of quadrature point
                        yk_hat = y_hat_NA(ui,z);

                        % Get x-component of displacement of node uk
                        vk = v(uk);
                        % Get y-component of displacement of node uk
                        wk = w(uk);

                        % Read per-bond quantities for bond ui-uk
                        IFk = IF_NA(ui,z);         % Influence function
                        Vk  = V_NA(ui,z);          % Neighbor area
                        Rk  = r_hat_NA(ui,z);      % Reference bond length

                        % Compute deformed bond components and length
                        rxk = xk_hat-xi + vk-vi;   % x-component of current relative position
                        ryk = yk_hat-yi + wk-wi;   % y-component of current relative position
                        rk  = sqrt(rxk^2 + ryk^2); % current bond length

                        % Compute bond stretch
                        sk = (rk - Rk)/Rk;

                        % ------------------------------------------------
                        %        Compute internal force density
                        % ------------------------------------------------

                        % Compute bond pairwise force magnitude
                        fk = c*IFk*sk;

                        % Compute bond pairwise force components
                        fvk = fk*rxk/rk;
                        fwk = fk*ryk/rk;

                        % Update internal force function components of node ui
                        Fv(ui) = Fv(ui) + fvk*Vk;
                        Fw(ui) = Fw(ui) + fwk*Vk;

                        % ------------------------------------------------
                        %      Compute macroelastic energy density
                        % ------------------------------------------------

                        % Compute bond pairwise potential
                        wk = 0.5*c*IFk*(sk^2)*Rk;

                        % Update macroelastic energy density of node ui
                        W(ui) = W(ui) + 0.5*wk*Vk;

                    end
                end
            end

        else

            error('flag_RDUG should be 0 or 1.')

        end

    % --------------------------------------------------------------------
    %                        Anisotropic model
    % --------------------------------------------------------------------
    elseif strcmp(model,'Anisotropic')
    
        % Number of nodes
        Nnodes = length(xx);
    
        % Initialize internal force density components and macroelastic energy density arrays
        Fv = zeros(Nnodes,1); % Array of x-components of internal force density for all nodes
        Fw = zeros(Nnodes,1); % Array of y-components of internal force density for all nodes
        W  = zeros(Nnodes,1); % Array of macroelastic energy density for all nodes
   
        % Find maximum number of neighbors any node can have
        zmax = length(u_NA(1,:));
    
        if flag_RDUG == 1

            % --------------------------------------------------------
            %              Computation for uniform grid
            %                 over rectangular domain
            % --------------------------------------------------------

            % Loop through all nodes ui
            for ui = 1:Nnodes

                % Get x-coordinate of node ui
                xi = xx(ui);
                % Get y-coordinate of node ui
                yi = yy(ui);

                % Get x-component of displacement of node ui
                vi = v(ui);
                % Get y-component of displacement of node ui
                wi = w(ui);

                % Loop over neighbors of node ui
                for z = 1:zmax

                    % Get neighbor cell number
                    uk = u_NA(ui,z);

                    % Only consider neighbor cells with a larger number
                    % than the source node to compute the bond pairwise force
                    % and pairwise potential once per bond
                    if uk > ui

                        % Get x-coordinate of quadrature point
                        xk_hat = x_hat_NA(ui,z);
                        % Get y-coordinate of quadrature point
                        yk_hat = y_hat_NA(ui,z);

                        % Get x-component of displacement of node uk
                        vk = v(uk);
                        % Get y-component of displacement of node uk
                        wk = w(uk);

                        % Read per-bond quantities for bond ui-uk
                        IFk = IF_NA(ui,z);         % Influence function
                        Vk  = V_NA(ui,z);          % Neighbor area
                        Rk  = r_hat_NA(ui,z);      % Reference bond length
                        Ck = c(ui,z);              % Micromodulus function

                        % Compute deformed bond components and length
                        rxk = xk_hat-xi + vk-vi;   % x-component of current relative position
                        ryk = yk_hat-yi + wk-wi;   % y-component of current relative position
                        rk  = sqrt(rxk^2 + ryk^2); % current bond length

                        % Compute bond stretch
                        sk = (rk - Rk)/Rk;

                        % ------------------------------------------------
                        %        Compute internal force density
                        % ------------------------------------------------

                        % Compute bond pairwise force magnitude
                        fk = Ck*IFk*sk;

                        % Compute bond pairwise force components
                        fvk = fk*rxk/rk;
                        fwk = fk*ryk/rk;

                        % Update internal force function components of node ui
                        Fv(ui) = Fv(ui) + fvk*Vk;
                        Fw(ui) = Fw(ui) + fwk*Vk;

                        % Update internal force function components of node uk
                        Fv(uk) = Fv(uk) - fvk*Vk;
                        Fw(uk) = Fw(uk) - fwk*Vk;

                        % ------------------------------------------------
                        %      Compute macroelastic energy density
                        % ------------------------------------------------

                        % Compute bond pairwise potential
                        wk = 0.5*Ck*IFk*(sk^2)*Rk;

                        % Update macroelastic energy density of node ui
                        W(ui) = W(ui) + 0.5*wk*Vk;

                        % Update macroelastic energy density of node uk
                        W(uk) = W(uk) + 0.5*wk*Vk;

                    end
                end
            end

        elseif flag_RDUG == 0    

            % --------------------------------------------------------
            %              Computation for general grids
            % --------------------------------------------------------

            % Loop through all nodes ui
            for ui = 1:Nnodes

                % Get x-coordinate of node ui
                xi = xx(ui);
                % Get y-coordinate of node ui
                yi = yy(ui);

                % Get x-component of displacement of node ui
                vi = v(ui);
                % Get y-component of displacement of node ui
                wi = w(ui);

                % Loop over neighbors of node ui
                for z = 1:zmax

                    % Get neighbor cell number
                    uk = u_NA(ui,z);

                    % Only consider neighbor cells
                    if uk > 0

                        % Get x-coordinate of quadrature point
                        xk_hat = x_hat_NA(ui,z);
                        % Get y-coordinate of quadrature point
                        yk_hat = y_hat_NA(ui,z);

                        % Get x-component of displacement of node uk
                        vk = v(uk);
                        % Get y-component of displacement of node uk
                        wk = w(uk);

                        % Read per-bond quantities for bond ui-uk
                        IFk = IF_NA(ui,z);         % Influence function
                        Vk  = V_NA(ui,z);          % Neighbor area
                        Rk  = r_hat_NA(ui,z);      % Reference bond length
                        Ck  = c(ui,z);             % Micromodulus function

                        % Compute deformed bond components and length
                        rxk = xk_hat-xi + vk-vi;   % x-component of current relative position
                        ryk = yk_hat-yi + wk-wi;   % y-component of current relative position
                        rk  = sqrt(rxk^2 + ryk^2); % current bond length

                        % Compute bond stretch
                        sk = (rk - Rk)/Rk;

                        % ------------------------------------------------
                        %        Compute internal force density
                        % ------------------------------------------------

                        % Compute bond pairwise force magnitude
                        fk = Ck*IFk*sk;

                        % Compute bond pairwise force components
                        fvk = fk*rxk/rk;
                        fwk = fk*ryk/rk;

                        % Update internal force function components of node ui
                        Fv(ui) = Fv(ui) + fvk*Vk;
                        Fw(ui) = Fw(ui) + fwk*Vk;

                        % ------------------------------------------------
                        %      Compute macroelastic energy density
                        % ------------------------------------------------

                        % Compute bond pairwise potential
                        wk = 0.5*Ck*IFk*(sk^2)*Rk;

                        % Update macroelastic energy density of node ui
                        W(ui) = W(ui) + 0.5*wk*Vk;

                    end
                end
            end

        else

            error('flag_RDUG should be 0 or 1.')

        end

    elseif strcmp(model, 'Classical')
        if flag_RDUG ~= 1
            error('Classical elasticity model requires rectangular domain with uniform grid.');
        end

        % Number of nodes
        Nnodes = length(xx);
        Ny = round(Nnodes / Nx);
    
        % Initialize internal force density components and macroelastic energy density arrays
        Fv = zeros(Nnodes,1); % Array of x-components of internal force density for all nodes
        Fw = zeros(Nnodes,1); % Array of y-components of internal force density for all nodes
        W  = zeros(Nnodes,1); % Array of macroelastic energy density for all nodes
        
        indices_x = zeros(Nnodes, 1);
        indices_y = zeros(Nnodes, 1);
        for i = 1:Ny
            for j = 1:Nx
                c = (i-1) * Nx + j;
                indices_y(c) = i;
                indices_x(c) = j;
            end
        end
        indices = sub2ind([Nx, Ny], indices_x, indices_y);
        
        dx = xx(2) - xx(1);
        dy = yy(1 + Nx) - yy(1);

        v_ = zeros(Nx, Ny);
        w_ = zeros(Nx, Ny);
        v_(indices) = v;
        w_(indices) = w;

        d2v_dx2 = zeros(Nx, Ny);
        d2v_dx2(2:end-1, :) = (v_(3:end, :) - 2*v_(2:end-1, :) + v_(1:end-2, :)) ./ (dx^2);
        d2v_dx2(1, :) = (v_(2, :) - 2*v_(1, :) + 0*v_(end, :)) ./ (dx^2);
        d2v_dx2(end, :) = (0*v_(1, :) - 2*v_(end, :) + v_(end-1, :)) ./ (dx^2);

        d2w_dy2 = zeros(Nx, Ny);
        d2w_dy2(:, 2:end-1) = (w_(:, 3:end) - 2*w_(:, 2:end-1) + w_(:, 1:end-2)) ./ (dy^2);
        d2w_dy2(:, 1) = (w_(:, 2) - 2*w_(:, 1) + 0*w_(:, end)) ./ (dy^2);
        d2w_dy2(:, end) = (0*w_(:, 1) - 2*w_(:, end) + w_(:, end-1)) ./ (dy^2);

        lap_v = zeros(Nx, Ny);
        lap_v(2:end-1, 2:end-1) = (v_(3:end, 3:end) + v_(1:end-2, 1:end-2) + v_(1:end-2, 3:end) + v_(3:end, 1:end-2) - 4*v_(2:end-1, 2:end-1)) / (2*dx^2);
        lap_v(1, 2:end-1) = (v_(2, 3:end) + 0*v_(end, 1:end-2) + 0*v_(end, 3:end) + v_(2, 1:end-2) - 4*v_(1, 2:end-1)) / (2*dx^2);
        lap_v(end, 2:end-1) = (0*v_(1, 3:end) + v_(end-1, 1:end-2) + v_(end-1, 3:end) + 0*v_(1, 1:end-2) - 4*v_(end, 2:end-1)) / (2*dx^2);
        lap_v(2:end-1, 1) = (v_(3:end, 2) + 0*v_(1:end-2, end) + v_(1:end-2, 2) + 0*v_(3:end, end) - 4*v_(2:end-1, 1)) / (2*dx^2);
        lap_v(2:end-1, end) = (0*v_(3:end, 1) + v_(1:end-2, end-1) + 0*v_(1:end-2, 1) + v_(3:end, end-1) - 4*v_(2:end-1, end)) / (2*dx^2);
        lap_v(1, 1) = (v_(2, 2) + 0*v_(end, end) + 0*v_(end, 2) + 0*v_(2, end) - 4*v_(1, 1)) / (2*dx^2);
        lap_v(1, end) = (0*v_(2, 1) + 0*v_(end, end-1) + 0*v_(end, 1) + v_(2, end-1) - 4*v_(1, end)) / (2*dx^2);
        lap_v(end, 1) = (0*v_(1, 2) + 0*v_(end-1, end) + v_(end-1, 2) + 0*v_(1, end) - 4*v_(end, 1)) / (2*dx^2);
        lap_v(end, end) = (0*v_(1, 1) + v_(end-1, end-1) + 0*v_(end-1, 1) + 0*v_(1, end-1) - 4*v_(end, end)) / (2*dx^2);

        lap_w = zeros(Nx, Ny);
        lap_w(2:end-1, 2:end-1) = (w_(3:end, 3:end) + w_(1:end-2, 1:end-2) + w_(1:end-2, 3:end) + w_(3:end, 1:end-2) - 4*w_(2:end-1, 2:end-1)) / (2*dx^2);
        lap_w(1, 2:end-1) = (w_(2, 3:end) + 0*w_(end, 1:end-2) + 0*w_(end, 3:end) + w_(2, 1:end-2) - 4*w_(1, 2:end-1)) / (2*dx^2);
        lap_w(end, 2:end-1) = (0*w_(1, 3:end) + w_(end-1, 1:end-2) + w_(end-1, 3:end) + 0*w_(1, 1:end-2) - 4*w_(end, 2:end-1)) / (2*dx^2);
        lap_w(2:end-1, 1) = (w_(3:end, 2) + 0*w_(1:end-2, end) + w_(1:end-2, 2) + 0*w_(3:end, end) - 4*w_(2:end-1, 1)) / (2*dx^2);
        lap_w(2:end-1, end) = (0*w_(3:end, 1) + w_(1:end-2, end-1) + 0*w_(1:end-2, 1) + w_(3:end, end-1) - 4*w_(2:end-1, end)) / (2*dx^2);
        lap_w(1, 1) = (w_(2, 2) + 0*w_(end, end) + 0*w_(end, 2) + 0*w_(2, end) - 4*w_(1, 1)) / (2*dx^2);
        lap_w(1, end) = (0*w_(2, 1) + 0*w_(end, end-1) + 0*w_(end, 1) + w_(2, end-1) - 4*w_(1, end)) / (2*dx^2);
        lap_w(end, 1) = (0*w_(1, 2) + 0*w_(end-1, end) + w_(end-1, 2) + 0*w_(1, end) - 4*w_(end, 1)) / (2*dx^2);
        lap_w(end, end) = (0*w_(1, 1) + w_(end-1, end-1) + 0*w_(end-1, 1) + 0*w_(1, end-1) - 4*w_(end, end)) / (2*dx^2);
        
        d2v_dxdy = zeros(Nx, Ny);
        d2v_dxdy(2:end-1, 2:end-1) = (v_(3:end, 3:end) - v_(3:end, 1:end-2) - v_(1:end-2, 3:end) + v_(1:end-2, 1:end-2)) / (4*dx*dy);
        d2v_dxdy(1, 2:end-1) = (v_(2, 3:end) - v_(2, 1:end-2) - 0*v_(end, 3:end) + 0*v_(end, 1:end-2)) / (4*dx*dy);
        d2v_dxdy(end, 2:end-1) = (0*v_(1, 3:end) - 0*v_(1, 1:end-2) - v_(end-1, 3:end) + v_(end-1, 1:end-2)) / (4*dx*dy);
        d2v_dxdy(2:end-1, 1) = (v_(3:end, 2) - 0*v_(3:end, end) - v_(1:end-2, 2) + 0*v_(1:end-2, end)) / (4*dx*dy);
        d2v_dxdy(2:end-1, end) = (0*v_(3:end, 1) - v_(3:end, end-1) - 0*v_(1:end-2, 1) + v_(1:end-2, end-1)) / (4*dx*dy);
        d2v_dxdy(1, 1) = (v_(2, 2) - 0*v_(2, end) - 0*v_(end, 2) + 0*v_(end, end)) / (4*dx*dy);
        d2v_dxdy(1, end) = (0*v_(2, 1) - v_(2, end-1) - 0*v_(end, 1) + 0*v_(end, end-1)) / (4*dx*dy);
        d2v_dxdy(end, 1) = (0*v_(1, 2) - 0*v_(1, end) - v_(end-1, 2) + 0*v_(end-1, end)) / (4*dx*dy);
        d2v_dxdy(end, end) = (0*v_(1, 1) - 0*v_(1, end-1) - 0*v_(end-1, 1) + v_(end-1, end-1)) / (4*dx*dy);

        d2w_dxdy = zeros(Nx, Ny);
        d2w_dxdy(2:end-1, 2:end-1) = (w_(3:end, 3:end) - w_(3:end, 1:end-2) - w_(1:end-2, 3:end) + w_(1:end-2, 1:end-2)) / (4*dx*dy);
        d2w_dxdy(1, 2:end-1) = (w_(2, 3:end) - w_(2, 1:end-2) - 0*w_(end, 3:end) + 0*w_(end, 1:end-2)) / (4*dx*dy);
        d2w_dxdy(end, 2:end-1) = (0*w_(1, 3:end) - 0*w_(1, 1:end-2) - w_(end-1, 3:end) + w_(end-1, 1:end-2)) / (4*dx*dy);
        d2w_dxdy(2:end-1, 1) = (w_(3:end, 2) - 0*w_(3:end, end) - w_(1:end-2, 2) + 0*w_(1:end-2, end)) / (4*dx*dy);
        d2w_dxdy(2:end-1, end) = (0*w_(3:end, 1) - w_(3:end, end-1) - 0*w_(1:end-2, 1) + w_(1:end-2, end-1)) / (4*dx*dy);
        d2w_dxdy(1, 1) = (w_(2, 2) - 0*w_(2, end) - 0*w_(end, 2) + 0*w_(end, end)) / (4*dx*dy);
        d2w_dxdy(1, end) = (0*w_(2, 1) - w_(2, end-1) - 0*w_(end, 1) + 0*w_(end, end-1)) / (4*dx*dy);
        d2w_dxdy(end, 1) = (0*w_(1, 2) - 0*w_(1, end) - w_(end-1, 2) + 0*w_(end-1, end)) / (4*dx*dy);
        d2w_dxdy(end, end) = (0*w_(1, 1) - 0*w_(1, end-1) - 0*w_(end-1, 1) + w_(end-1, end-1)) / (4*dx*dy);

        Fv_ = (3/4*E) * (d2v_dx2 + d2w_dxdy + (1/2) * lap_v);
        Fw_ = (3/4*E) * (d2w_dy2 + d2v_dxdy + (1/2) * lap_w);
        Fv(:) = Fv_(indices);
        Fw(:) = Fw_(indices);

    elseif strcmp(model, 'Classical-periodic')
        if flag_RDUG ~= 1
            error('Classical elasticity model requires rectangular domain with uniform grid.');
        end

        % Number of nodes
        Nnodes = length(xx);
        Ny = round(Nnodes / Nx);
    
        % Initialize internal force density components and macroelastic energy density arrays
        Fv = zeros(Nnodes,1); % Array of x-components of internal force density for all nodes
        Fw = zeros(Nnodes,1); % Array of y-components of internal force density for all nodes
        W  = zeros(Nnodes,1); % Array of macroelastic energy density for all nodes
        
        indices_x = zeros(Nnodes, 1);
        indices_y = zeros(Nnodes, 1);
        for i = 1:Ny
            for j = 1:Nx
                c = (i-1) * Nx + j;
                indices_y(c) = i;
                indices_x(c) = j;
            end
        end
        indices = sub2ind([Nx, Ny], indices_x, indices_y);
        
        dx = xx(2) - xx(1);
        dy = yy(1 + Nx) - yy(1);

        v_ = zeros(Nx, Ny);
        w_ = zeros(Nx, Ny);
        v_(indices) = v;
        w_(indices) = w;

        d2v_dx2 = zeros(Nx, Ny);
        d2v_dx2(2:end-1, :) = (v_(3:end, :) - 2*v_(2:end-1, :) + v_(1:end-2, :)) ./ (dx^2);
        d2v_dx2(1, :) = (v_(2, :) - 2*v_(1, :) + v_(end, :)) ./ (dx^2);
        d2v_dx2(end, :) = (v_(1, :) - 2*v_(end, :) + v_(end-1, :)) ./ (dx^2);

        d2w_dy2 = zeros(Nx, Ny);
        d2w_dy2(:, 2:end-1) = (w_(:, 3:end) - 2*w_(:, 2:end-1) + w_(:, 1:end-2)) ./ (dy^2);
        d2w_dy2(:, 1) = (w_(:, 2) - 2*w_(:, 1) + w_(:, end)) ./ (dy^2);
        d2w_dy2(:, end) = (w_(:, 1) - 2*w_(:, end) + w_(:, end-1)) ./ (dy^2);

        lap_v = zeros(Nx, Ny);
        lap_v(2:end-1, 2:end-1) = (v_(3:end, 3:end) + v_(1:end-2, 1:end-2) + v_(1:end-2, 3:end) + v_(3:end, 1:end-2) - 4*v_(2:end-1, 2:end-1)) / (2*dx^2);
        lap_v(1, 2:end-1) = (v_(2, 3:end) + v_(end, 1:end-2) + v_(end, 3:end) + v_(2, 1:end-2) - 4*v_(1, 2:end-1)) / (2*dx^2);
        lap_v(end, 2:end-1) = (v_(1, 3:end) + v_(end-1, 1:end-2) + v_(end-1, 3:end) + v_(1, 1:end-2) - 4*v_(end, 2:end-1)) / (2*dx^2);
        lap_v(2:end-1, 1) = (v_(3:end, 2) + v_(1:end-2, end) + v_(1:end-2, 2) + v_(3:end, end) - 4*v_(2:end-1, 1)) / (2*dx^2);
        lap_v(2:end-1, end) = (v_(3:end, 1) + v_(1:end-2, end-1) + v_(1:end-2, 1) + v_(3:end, end-1) - 4*v_(2:end-1, end)) / (2*dx^2);
        lap_v(1, 1) = (v_(2, 2) + v_(end, end) + v_(end, 2) + v_(2, end) - 4*v_(1, 1)) / (2*dx^2);
        lap_v(1, end) = (v_(2, 1) + v_(end, end-1) + v_(end, 1) + v_(2, end-1) - 4*v_(1, end)) / (2*dx^2);
        lap_v(end, 1) = (v_(1, 2) + v_(end-1, end) + v_(end-1, 2) + v_(1, end) - 4*v_(end, 1)) / (2*dx^2);
        lap_v(end, end) = (v_(1, 1) + v_(end-1, end-1) + v_(end-1, 1) + v_(1, end-1) - 4*v_(end, end)) / (2*dx^2);

        lap_w = zeros(Nx, Ny);
        lap_w(2:end-1, 2:end-1) = (w_(3:end, 3:end) + w_(1:end-2, 1:end-2) + w_(1:end-2, 3:end) + w_(3:end, 1:end-2) - 4*w_(2:end-1, 2:end-1)) / (2*dx^2);
        lap_w(1, 2:end-1) = (w_(2, 3:end) + w_(end, 1:end-2) + w_(end, 3:end) + w_(2, 1:end-2) - 4*w_(1, 2:end-1)) / (2*dx^2);
        lap_w(end, 2:end-1) = (w_(1, 3:end) + w_(end-1, 1:end-2) + w_(end-1, 3:end) + w_(1, 1:end-2) - 4*w_(end, 2:end-1)) / (2*dx^2);
        lap_w(2:end-1, 1) = (w_(3:end, 2) + w_(1:end-2, end) + w_(1:end-2, 2) + w_(3:end, end) - 4*w_(2:end-1, 1)) / (2*dx^2);
        lap_w(2:end-1, end) = (w_(3:end, 1) + w_(1:end-2, end-1) + w_(1:end-2, 1) + w_(3:end, end-1) - 4*w_(2:end-1, end)) / (2*dx^2);
        lap_w(1, 1) = (w_(2, 2) + w_(end, end) + w_(end, 2) + w_(2, end) - 4*w_(1, 1)) / (2*dx^2);
        lap_w(1, end) = (w_(2, 1) + w_(end, end-1) + w_(end, 1) + w_(2, end-1) - 4*w_(1, end)) / (2*dx^2);
        lap_w(end, 1) = (w_(1, 2) + w_(end-1, end) + w_(end-1, 2) + w_(1, end) - 4*w_(end, 1)) / (2*dx^2);
        lap_w(end, end) = (w_(1, 1) + w_(end-1, end-1) + w_(end-1, 1) + w_(1, end-1) - 4*w_(end, end)) / (2*dx^2);
        
        d2v_dxdy = zeros(Nx, Ny);
        d2v_dxdy(2:end-1, 2:end-1) = (v_(3:end, 3:end) - v_(3:end, 1:end-2) - v_(1:end-2, 3:end) + v_(1:end-2, 1:end-2)) / (4*dx*dy);
        d2v_dxdy(1, 2:end-1) = (v_(2, 3:end) - v_(2, 1:end-2) - v_(end, 3:end) + v_(end, 1:end-2)) / (4*dx*dy);
        d2v_dxdy(end, 2:end-1) = (v_(1, 3:end) - v_(1, 1:end-2) - v_(end-1, 3:end) + v_(end-1, 1:end-2)) / (4*dx*dy);
        d2v_dxdy(2:end-1, 1) = (v_(3:end, 2) - v_(3:end, end) - v_(1:end-2, 2) + v_(1:end-2, end)) / (4*dx*dy);
        d2v_dxdy(2:end-1, end) = (v_(3:end, 1) - v_(3:end, end-1) - v_(1:end-2, 1) + v_(1:end-2, end-1)) / (4*dx*dy);
        d2v_dxdy(1, 1) = (v_(2, 2) - v_(2, end) - v_(end, 2) + v_(end, end)) / (4*dx*dy);
        d2v_dxdy(1, end) = (v_(2, 1) - v_(2, end-1) - v_(end, 1) + v_(end, end-1)) / (4*dx*dy);
        d2v_dxdy(end, 1) = (v_(1, 2) - v_(1, end) - v_(end-1, 2) + v_(end-1, end)) / (4*dx*dy);
        d2v_dxdy(end, end) = (v_(1, 1) - v_(1, end-1) - v_(end-1, 1) + v_(end-1, end-1)) / (4*dx*dy);

        d2w_dxdy = zeros(Nx, Ny);
        d2w_dxdy(2:end-1, 2:end-1) = (w_(3:end, 3:end) - w_(3:end, 1:end-2) - w_(1:end-2, 3:end) + w_(1:end-2, 1:end-2)) / (4*dx*dy);
        d2w_dxdy(1, 2:end-1) = (w_(2, 3:end) - w_(2, 1:end-2) - w_(end, 3:end) + w_(end, 1:end-2)) / (4*dx*dy);
        d2w_dxdy(end, 2:end-1) = (w_(1, 3:end) - w_(1, 1:end-2) - w_(end-1, 3:end) + w_(end-1, 1:end-2)) / (4*dx*dy);
        d2w_dxdy(2:end-1, 1) = (w_(3:end, 2) - w_(3:end, end) - w_(1:end-2, 2) + w_(1:end-2, end)) / (4*dx*dy);
        d2w_dxdy(2:end-1, end) = (w_(3:end, 1) - w_(3:end, end-1) - w_(1:end-2, 1) + w_(1:end-2, end-1)) / (4*dx*dy);
        d2w_dxdy(1, 1) = (w_(2, 2) - w_(2, end) - w_(end, 2) + w_(end, end)) / (4*dx*dy);
        d2w_dxdy(1, end) = (w_(2, 1) - w_(2, end-1) - w_(end, 1) + w_(end, end-1)) / (4*dx*dy);
        d2w_dxdy(end, 1) = (w_(1, 2) - w_(1, end) - w_(end-1, 2) + w_(end-1, end)) / (4*dx*dy);
        d2w_dxdy(end, end) = (w_(1, 1) - w_(1, end-1) - w_(end-1, 1) + w_(end-1, end-1)) / (4*dx*dy);

        Fv_ = (3/4*E) * (d2v_dx2 + d2w_dxdy + (1/2) * lap_v);
        Fw_ = (3/4*E) * (d2w_dy2 + d2v_dxdy + (1/2) * lap_w);
        Fv(:) = Fv_(indices);
        Fw(:) = Fw_(indices);
    else
        error('Invalid model.')

    end

end
