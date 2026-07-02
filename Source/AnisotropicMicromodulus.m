
% ========================================================================
% The function AnisotropicMicromodulus computes the anisotropic
% micromodulus function for a given material parameterized by its
% peridynamic tensor
% ========================================================================

% Input
% -----
% xi         : reference bond, i.e., x' - x
% m          : micromodulus normalization
% Lambda     : peridynamic tensor
%

% Output
% ------
% c          : value of micromodulus function c(xi)
%

function c = AnisotropicMicromodulus(xi, m, Lambda)
    xiMag = norm(xi);
    a = xiMag / m;
    xi = xi / xiMag;

    s = Lambda(1, 1, 1, 1) * xi(1)^4;
    s = s + 4 * Lambda(1, 1, 1, 2) * xi(1)^3 * xi(2);
    s = s + 6 * Lambda(1, 1, 2, 2) * xi(1)^2 * xi(2)^2;
    s = s + 4 * Lambda(2, 2, 1, 2) * xi(1) * xi(2)^3;
    s = s + Lambda(2, 2, 2, 2) * xi(2)^4;
    
    c = a * s;
end