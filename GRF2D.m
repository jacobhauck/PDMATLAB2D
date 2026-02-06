function f = GRF2D(num_modes, amplitude_fn)
    % Generates a zero-mean Gaussian random field defined on [0,1]^2
    % 
    % Parameters
    % ----------
    % num_modes: Number of Fourier modes (in each direction, + and -) used
    %            to simulate the field
    % amplitude_fn: Function with signature (i, j) -> (a), where i and j
    %               are the Fourier mode index arrays of arbitrary shape,
    %               and a is an array of amplitudes with the same shape as
    %               i and j
    %
    % Outputs
    % -------
    % f: Callable with signature (x, y) -> (f), where x and y are arrays
    %    of x and y coordinates (of any shape), and f is an array of the
    %    values of the field at x and y (with the same shape)
 
    [gx, gy] = meshgrid(-num_modes:num_modes, -num_modes:num_modes);
    % (2*num_modes + 1, 2*num_modes + 1) each
    % total_modes = (2*num_modes + 1)^2
    gx = reshape(gx, numel(gx), 1);  % (total_modes, 1)
    gy = reshape(gy, numel(gy), 1);  % (total_modes, 1)
    kx = (2*pi) * reshape(gx, numel(gx), 1);  % (total_modes, 1)
    ky = (2*pi) * reshape(gy, numel(gy), 1);  % (total_modes, 1)

    a = amplitude_fn(gx, gy);  % (total_modes, 1)
    sinCoef = a .* randn(size(a));  % (total_modes, 1)
    cosCoef = a .* randn(size(a));  % (total_modes, 1)
    
    f = @(x, y) evalFourierSeries(x, y, kx, ky, sinCoef, cosCoef);
end

function fVal = evalFourierSeries(x, y, kx, ky, sinCoef, cosCoef)
    % Evaluates a Fourier series at given points
    %
    % Parameters
    % ----------
    % x: (*shape) array of x values
    % y: (*shape) array of y values
    % kx: (total_modes, 1) array of wave vector x components
    % ky: (total_modes, 1) array of wave vector y components
    % sinCoef: (total_modes, 1) array of coefficients of sine terms
    % cosCoef: (total_modes, 1) array of coefficients of cosine terms
    %
    % Outputs
    % -------
    % fVal: (*shape) array of of values of fourier series at x and y
    
    shape = size(x);
    x = reshape(x, numel(x), 1);  % (n, 1)
    y = reshape(y, numel(y), 1);  % (n, 1)

    phase = x * kx' + y * ky';  % (n, total_modes)
    fVal = sin(phase) * sinCoef + cos(phase) * cosCoef;  % (n, 1)
    fVal = reshape(fVal, shape);
end
