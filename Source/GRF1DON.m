classdef GRF1DON
    properties
        % ===== Public properties =====
        % Number of modes (cos and sin)
        numModes = 8

        % Amplitude decay
        amplitudeFn = @(g) ones(size(g))
        
        % Domain left
        x0 = 0

        % Domain right
        x1 = 1

        % Optional custom stream of random numbers
        rngStream = []

        % ===== Private properties =====
        ak
        an
        aAmplitude
        aCov

        bk
        bn
        bAmplitude
        bCov

        cov
    end
    
    methods
        function self = GRF1DON(numModes, amplitudeFn, domain, rngStream)
            % A zero-mean Gaussian random field defined on [x0, x1] using
            % an orthonormal basis
            % 
            % Parameters
            % ----------
            % numModes: Number of Fourier modes used to simulate the field
            % amplitudeFn: Function with signature (g) -> (a), where g
            %              contains the integer wave number of the mode
            %              and a is an array of amplitudes with the same 
            %              shape as g
            % domain: (2, 1) column vector [x0; x1] defining the
            %         domain [x0, x1] of the GRF
            % rngStream: Optional random number stream to use to generate
            %            coefficients

            self.numModes = numModes;
            self.amplitudeFn = amplitudeFn;
            self.x0 = domain(1);
            self.x1 = domain(2);
            
            if nargin == 4
                self.rngStream = rngStream;
            end

            gCos = (0:numModes)';  % (numModes + 1, 1)
            gSin = (1:numModes)';  % (numModes, 1)
            
            scale = sqrt(self.x1 - self.x0);

            amplitudeCos = amplitudeFn(gCos);  % (numModes + 1, 1)
            nCos = [scale; scale * ones(size(gSin)) / sqrt(2)];  % (numModes + 1, 1)
            amplitudeSin = amplitudeFn(gSin);  % (numModes, 1)
            nSin = scale * ones(size(gSin)) / sqrt(2);  % (numModes, 1)
            
            % all shape (numModes + 1, 1)
            self.ak = (2*pi) * gCos;
            self.an = nCos;
            self.aAmplitude = amplitudeCos;
            self.aCov = self.aAmplitude .^ 2;
        
            % all shape (numModes, 1)
            self.bk = (2*pi) * gSin;
            self.bn = nSin;
            self.bAmplitude = amplitudeSin;
            self.bCov = self.bAmplitude .^ 2;

            self.cov = zeros(self.basisDimension(), self.basisDimension());
            self.cov(1:length(self.ak), 1:length(self.ak)) = diag(self.aCov);
            start = length(self.ak) + 1;
            self.cov(start:end, start:end) = diag(self.bCov);
        end

        function f = generate(self, rngStream)
            % Generates a sample from the GRF as a function
            % 
            % Parameters
            % ----------
            % rngStream: Optional stream of random numbers to use instead
            %            of the default or self.rngStream
            %
            % Return
            % ------
            % f: function (x) -> f(x) sampled from the GRF. x may
            % have any shape, and f(x) will have the same shape, with f
            % applied elementwise.
            
            if nargin == 1
                a = randn(rngStream, size(self.ak));  % (numModes + 1, 1)
                b = randn(rngStream, size(self.bk));  % (numModes, 1)
            elseif isempty(self.rngStream)
                a = randn(size(self.ak));  % (numModes + 1, 1)
                b = randn(size(self.bk));  % (numModes, 1)
            else
                a = randn(self.rngStream, size(self.ak));  % (numModes + 1, 1)
                b = randn(self.rngStream, size(self.bk));  % (numModes, 1)
            end
            
            coef = [
                a .* self.aAmplitude;
                b .* self.bAmplitude;
            ];  % (n, 1)

            f = @(x) self.eval(x, coef);
        end

        function fVal = eval(self, x, coef)
            % Evaluates the GRF Fourier series for particular coefficients
            %
            % Parameters
            % ----------
            % x: (shape) x coordinates
            % coef: (n, 1) vector of coefficients (in order a, b, c, d)
            %
            % Return
            % ------
            % fVal: (shape) Fourier series evaluated (elementwise) at 
            %       (x)
            a = coef(1 : length(self.ak));  % (numModes + 1, 1)
            start = length(self.ak) + 1;
            b = coef(start:end);  % (numModes, 1)

            aCoef = a ./ self.an;  % (numModes + 1, 1)
            bCoef = b ./ self.bn;  % (numModes, 1)

            originalSize = size(x);
            x = reshape(x, 1, []);  % (1, numPoints)
            
            % Transform to standard domain
            x = (x - self.x0) ./ (self.x1 - self.x0);  % (1, numPoints)

            a = aCoef' * cos(self.ak * x);  % (1, numPoints)
            b = bCoef' * sin(self.bk * x);  % (1, numPoints)
        
            fVal = reshape(a + b, originalSize);  % (*shape)
        end

        function fVal = evalBasis(self, x, i)
            % Evaluates a single GRF Fourier series basis function
            %
            % Parameters
            % ----------
            % x: (shape) x coordinates
            % i: index of basis function to evaluate
            %
            % Return
            % ------
            % fVal: (shape) basis function i evaluated (elementwise) at 
            %       (x, y)

            originalSize = size(x);
            x = reshape(x, 1, []);  % (1, numPoints)
            
            % Transform to standard domain
            x = (x - self.x0) ./ (self.x1 - self.x0);  % (1, numPoints)

            if i <= length(self.ak)
                fVal = cos(self.ak(i) * x) ./ self.an(i);
            else
                i = i - length(self.ak);
                fVal = sin(self.bk(i) * x) ./ self.bn(i);
            end

            fVal = reshape(fVal, originalSize);  % (*shape)
        end

        function showBasis(self)
            % Shows the basis functions of the GRF
            x = linspace(self.x0, self.x1, 128)';
            
            counter = 1;
            total = self.basisDimension();
            
            % Show a basis functions
            for i = 1:length(self.ak)
                f = self.evalBasis(x, counter);
                plot(x, f);
                g = round(self.ak(i) / (2*pi));
                title(sprintf("a (cos) basis func %d, g = %d", i, g));
                input(sprintf("Press Enter to view next (%d/%d)", counter, total));
                xlabel('x');
                counter = counter + 1;
            end

            % Show b basis functions
            for i = 1:length(self.bk)
                f = self.evalBasis(x, counter);
                plot(x, f);
                g = round(self.bk(i) / (2*pi));
                title(sprintf("b (sin) basis func %d, g = %d", i, g));
                xlabel('x');
                input(sprintf("Press Enter to view next (%d/%d)", counter, total));
                counter = counter + 1;
            end
        end

        function validateOrthonormality(self)
            % Validates that the basis functions are orthonormal
            x = linspace(self.x0, self.x1, 129)';
            
            % Construct basis functions
            basis = zeros(self.basisDimension(), size(x, 1));
            for i = 1:self.basisDimension()
                basis(i, 1:end) = self.evalBasis(x, i);
            end

            % Evaluate inner products
            gramian = zeros(self.basisDimension(), self.basisDimension());
            for i = 1:self.basisDimension()
                for j = 1:self.basisDimension()
                    basisi = reshape(basis(i, 1:end), [], 1);
                    basisj = reshape(basis(j, 1:end), [], 1);
                    dx = (self.x1 - self.x0) / 128;
                    gramian(i, j) = (basisi' * basisj) * dx;
                end
            end

            fprintf("gramian\n");
            gramian

            fprintf( ...
                "Max deviation from I: %f\n", ...
                max(max(abs(gramian - eye(self.basisDimension())))) ...
            );
        end

        function uCoef = coefficients(self, u)
            % Computes the coefficients of a given function in the GRF's
            % basis
            %
            % Parameters
            % ----------
            % u: function (x) -> u(x) that applies the desired
            %    function elementwise to array x
            %
            % Return
            % ------
            % uCoef: (n, 1) coefficients of the given function in the GRF
            %        basis

            x = linspace(self.x0, self.x1, 129)';  % (numPoints, 1)
            uVal = u(x);  % (numPoints, 1)

            uCoef = zeros(self.basisDimension(), 1);
            for i = 1:self.basisDimension()
                basisi = self.evalBasis(x, y, i);  % (numPoints, 1)
                dx = (self.x1 - self.x0) / 128;
                uCoef(i) = (basisi' * uVal) * dx;
            end
        end

        function n = basisDimension(self)
            % Gets the dimension of the basis used by the GRF
            %
            % Return
            % ------
            % n: Dimension of the basis

            n = length(self.ak) + length(self.bk);
        end
    end
end