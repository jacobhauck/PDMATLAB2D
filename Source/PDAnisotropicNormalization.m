
% ========================================================================
% The function PDAnisotropicNormalization computes the micromodulus
% function normaliization factor for the Anisotropic PMB model.
% ========================================================================

% Input
% -----
% omega       : influence function order indicator
% del         : horizon

% Output
% ------
% m           : micromodulus normalization factor

function m = PDAnisotropicNormalization(omega, del)
    switch omega
        case 0
            m = pi * del^4 / 2;
        case 0.5
            m = pi * del^4 / 2;
        case 1
            m = pi * del^4 / 10;
        case 3
            m = pi * del^4 / 14;
        case 5
            m = 5 * pi * del^4 / 84;
        case 7
            m = 7 * pi * del^4 / 132;
        otherwise
            error("Invalid omega");
    end
end