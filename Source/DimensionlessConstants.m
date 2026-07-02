function [E, Go, del] = DimensionlessConstants(omega, n, so, model, PlanarModel)
    del = 1/n;  % Target delta

    % --------------------------------------------------------------------
    %                           GPMB model
    % --------------------------------------------------------------------

    if strcmp(model,'GPMB')
    
        % ----------------------------------------------------------------
        %                         Plane strain
        % ----------------------------------------------------------------

        if strcmp(PlanarModel,'PlaneStrain')
 
              % Material constants for each influence function (IF)
              switch omega
                 case 0
                     % Constant IF
                     E = 5 / 48 / n;
                     Go = so^2 * 12 * E * del / (5 * pi);
                 case 0.5
                     % Piecewise constant IF
                     E = 5 / 48 / n;
                     Go = so^2 * 12 * E * del / (5 * pi);
                 case 1
                     % Piecewise linear IF
                     E = 5 / 192 / n;
                     Go = so^2 * 48 * E * del / (25 * pi);
                 case 3
                     % Piecewise cubic IF
                     E = 1 / 48 / n;
                     Go = so^2 * 12 * E * del / (7 * pi);
                 case 5
                     % Piecewise quintic IF
                     E = 25 / 1344 / n;
                     Go = so^2 * 8 * E * del / (5 * pi);
                 case 7
                     % Piecewise septic IF
                     E = 5 / 288 / n;
                     Go = so^2 * 84 * E * del / (55 * pi);
                 otherwise
                     error('Invalid omega.');
              end

        % ----------------------------------------------------------------
        %                         Plane stress
        % ----------------------------------------------------------------
    
        elseif strcmp(PlanarModel,'PlaneStress')
    
            % Material constants for each influence function (IF)
            switch omega
                case 0
                    % Constant IF
                    E = 1 / 9 / n;
                    Go = so^2 * 9 * E * del / (4 * pi);
                case 0.5
                    % Piecewise constant IF
                    E = 1 / 9 / n;
                    Go = so^2 * 9 * E * del / (4 * pi);
                case 1
                    % Piecewise linear IF
                    E = 1 / 36 / n;
                    Go = so^2 * 9 * E * del / (5 * pi);
                case 3
                    % Piecewise cubic IF
                    E = 1 / 45 / n;
                    Go = so^2 * 45 * E * del / (28 * pi);
                case 5
                    % Piecewise quintic IF
                    E = 5 / 252 / n;
                    Go = so^2 * 3 * E * del / (2 * pi);
                case 7
                    % Piecewise septic IF
                    E = 1 / 54 / n;
                    Go = so^2 * 63 * E * del / (44 * pi);
                otherwise
                    error('Invalid omega.');
            end
    
        else
    
            error('Invalid PlanarModel.')
    
        end
    
    else
    
        error('Invalid model.')
    
    end

end