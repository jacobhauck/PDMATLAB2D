function f = GRF1D(num_modes, amplitude_fn)
    % Generates a zero-mean Gaussian random field defined on [0,1]
    % 
    % Parameters
    % ----------
    % num_modes: Number of Fourier modes (in each direction, + and -) used
    %            to simulate the field
    % amplitude_fn: Function with signature (i) -> (a), where i is the
    %               Fourier mode index array of arbitrary shape,
    %               and a is an array of amplitudes with the same shape as
    %               i
    %
    % Outputs
    % -------
    % f: Callable with signature (x) -> (f), where x is an array
    %    of x coordinates (of any shape), and f is an array of the
    %    values of the field at x (with the same shape)
 
    gx = (-num_modes:num_modes)';
    % (2*num_modes + 1, 1) 
    kx = (2*pi) * gx;  % (total_modes, 1)

    a = amplitude_fn(gx);  % (total_modes, 1)
    sinCoef = a .* randn(size(a));  % (total_modes, 1)
    cosCoef = a .* randn(size(a));  % (total_modes, 1)
    
    f = @(x) evalFourierSeries(x, kx, sinCoef, cosCoef);
end

function fVal = evalFourierSeries(x, kx, sinCoef, cosCoef)
    % Evaluates a Fourier series at given points
    %
    % Parameters
    % ----------
    % x: (*shape) array of x values
    % kx: (total_modes, 1) array of wave vector x components
    % sinCoef: (total_modes, 1) array of coefficients of sine terms
    % cosCoef: (total_modes, 1) array of coefficients of cosine terms
    %
    % Outputs
    % -------
    % fVal: (*shape) array of of values of Fourier series at x
    
    shape = size(x);
    x = reshape(x, numel(x), 1);  % (n, 1)

    phase = x * kx';  % (n, total_modes)
    fVal = sin(phase) * sinCoef + cos(phase) * cosCoef;  % (n, 1)
    fVal = reshape(fVal, shape);
end
