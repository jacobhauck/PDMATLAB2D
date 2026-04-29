classdef GRF2DON
    properties
        % ===== Public properties =====
        % Number of modes in each direction
        numModes = 8

        % Amplitude decay
        amplitudeFn = @(gx, gy) ones(size(gx))
        
        % Domain left
        x0 = 0

        % Domain right
        x1 = 1

        % Domain bottom
        y0 = 0

        % Domain top
        y1 = 1

        % Optional custom stream of random numbers
        rngStream = []

        % ===== Private properties =====
        akx
        aky
        an
        aAmplitude
        aCov

        bkx
        bky
        bn
        bAmplitude
        bCov
        
        ckx
        cky
        cn
        cAmplitude
        cCov
        
        dkx
        dky
        dn
        dAmplitude
        dCov

        cov
    end

    methods
        function self = GRF2DON(numModes, amplitudeFn, domain, rngStream)
            % A zero-mean Gaussian random field defined on 
            % [x0, x1] x [y0, y1] using an orthonormal Fourier basis
            % 
            % Parameters
            % ----------
            % numModes: Number of Fourier modes (in each direction) used
            %            to simulate the field
            % amplitudeFn: Function with signature (gx, gy) -> (a), where gx,
            %               and gy are the integer wave numbers of the mode
            %               and a is an array of amplitudes with the same shape as
            %               gx, gy
            % domain: (4, 1) column vector [x0; y0; x1; y1] defining the
            %         domain [x0, x1] x [y0, y1] of the GRF
            % rngStream: Optional random number stream to use to generate
            %            coefficients

            self.numModes = numModes;
            self.amplitudeFn = amplitudeFn;
            self.x0 = domain(1);
            self.x1 = domain(3);
            self.y0 = domain(2);
            self.y1 = domain(4);
            
            if nargin == 4
                self.rngStream = rngStream;
            end

            gCos = (0:numModes)';  % (numModes + 1, 1)
            gSin = (1:numModes)';  % (numModes, 1)
            
            [gx, gy] = meshgrid(gSin, gSin);  % (numModes, numModes) each
            gxAll = reshape(gx, numModes^2, 1);  % (numModes^2, 1)
            gyAll = reshape(gy, numModes^2, 1);  % (numModes^2, 1)
            
            scale = sqrt((self.x1 - self.x0) * (self.y1 - self.y0));

            amplitudeAll = amplitudeFn(gxAll, gyAll);  % (numModes^2, 1)
            nAll = scale * ones(size(amplitudeAll)) / 2;  % (numModes^2, 1)
            amplitudeCos = amplitudeFn(gCos, zeros(size(gCos)));  % (numModes + 1, 1)
            nCos = [scale; scale * ones(size(gSin)) / sqrt(2)];  % (numModes + 1, 1)
            amplitudeSin = amplitudeFn(gSin, zeros(size(gSin)));  % (numModes, 1)
            nSin = scale * ones(size(gSin)) / sqrt(2);  % (numModes, 1)
            
            % all shape ((numModes + 1)^2, 1)
            self.akx = (2*pi) * [gxAll; gCos; zeros(size(gSin))];  
            self.aky = (2*pi) * [gyAll; zeros(size(gCos)); gSin];
            self.an = [nAll; nCos; nSin];
            self.aAmplitude = [amplitudeAll; amplitudeCos; amplitudeSin];
            self.aCov = self.aAmplitude .^ 2;
        
            % all shape (numModes^2 + numModes, 1)
            self.bkx = (2*pi) * [gxAll; zeros(size(gSin))];  
            self.bky = (2*pi) * [gyAll; gSin];
            self.bn = [nAll; nSin];
            self.bAmplitude = [amplitudeAll; amplitudeSin];
            self.bCov = self.bAmplitude .^ 2;

            % all shape (numModes^2 + numModes, 1)
            self.ckx = (2*pi) * [gxAll; gSin];
            self.cky = (2*pi) * [gyAll; zeros(size(gSin))];
            self.cn = [nAll; nSin];
            self.cAmplitude = [amplitudeAll; amplitudeSin];
            self.cCov = self.cAmplitude .^ 2;

            % all shape (numModes^2, 1)
            self.dkx = (2*pi) * gxAll;
            self.dky = (2*pi) * gyAll;
            self.dn = nAll;
            self.dAmplitude = amplitudeAll;
            self.dCov = self.dAmplitude .^ 2;

            self.cov = zeros(self.basisDimension(), self.basisDimension());
            self.cov(1:length(self.akx), 1:length(self.akx)) = diag(self.aCov);
            start = length(self.akx) + 1;
            self.cov(start : start + length(self.bkx) - 1, start : start + length(self.bkx) - 1) = diag(self.bCov);
            start = start + length(self.bkx);
            self.cov(start : start + length(self.ckx) - 1, start : start + length(self.ckx) - 1) = diag(self.cCov);
            start = start + length(self.ckx);
            self.cov(start:end, start:end) = diag(self.dCov);
        end

        function f = generate(self, rngStream)
            % Generates a sample from the GRF as a function
            % 
            % Parameters
            % ----------
            % rngStream: Optional stream of random numbers to use instead
            %            of default or self.rngStream
            %
            % Return
            % ------
            % f: function (x, y) -> f(x, y) sampled from the GRF. x, y may
            % have any shape, and f(x, y) will have the same shape, with f
            % applied elementwise.
            
            if nargin == 1
                a = randn(rngStream, size(self.akx));  % ((numModes + 1)^2, 1)
                b = randn(rngStream, size(self.bkx));  % (numModes^2 + numModes, 1)
                c = randn(rngStream, size(self.ckx));  % (numModes^2 + numModes, 1)
                d = randn(rngStream, size(self.dkx));  % (numModes^2, 1)
            elseif isempty(self.rngStream)
                a = randn(size(self.akx));  % ((numModes + 1)^2, 1)
                b = randn(size(self.bkx));  % (numModes^2 + numModes, 1)
                c = randn(size(self.ckx));  % (numModes^2 + numModes, 1)
                d = randn(size(self.dkx));  % (numModes^2, 1)
            else
                a = randn(self.rngStream, size(self.akx));  % ((numModes + 1)^2, 1)
                b = randn(self.rngStream, size(self.bkx));  % (numModes^2 + numModes, 1)
                c = randn(self.rngStream, size(self.ckx));  % (numModes^2 + numModes, 1)
                d = randn(self.rngStream, size(self.dkx));  % (numModes^2, 1)
            end
            
            coef = [
                a .* self.aAmplitude;
                b .* self.bAmplitude;
                c .* self.cAmplitude;
                d .* self.dAmplitude
            ];  % (n, 1)

            f = @(x, y) self.eval(x, y, coef);
        end

        function fVal = eval(self, x, y, coef)
            % Evaluates the GRF Fourier series for particular coefficients
            %
            % Parameters
            % ----------
            % x: (shape) x coordinates
            % y: (shape) y coordinates
            % coef: (n, 1) vector of coefficients (in order a, b, c, d)
            %
            % Return
            % ------
            % fVal: (shape) Fourier series evaluated (elementwise) at 
            %       (x, y)

            start = 1;
            a = coef(start : length(self.akx));  % ((numModes + 1)^2, 1)
            start = start + length(a);
            b = coef(start : start + length(self.bkx) - 1);  % ((numModes^2 + numModes, 1)
            start = start + length(b);
            c = coef(start : start + length(self.ckx) - 1);  % ((numModes^2 + numModes, 1)
            start = start + length(c);
            d = coef(start:end);  % ((numModes^2, 1)

            aCoef = a ./ self.an;  % ((numModes + 1)^2, 1)
            bCoef = b ./ self.bn;  % ((numModes^2 + numModes, 1)
            cCoef = c ./ self.cn;  % ((numModes^2 + numModes, 1)
            dCoef = d ./ self.dn;  % ((numModes^2, 1)

            originalSize = size(x);
            x = reshape(x, 1, []);  % (1, numPoints)
            y = reshape(y, 1, []);  % (1, numPoints)
            
            % Transform to standard domain
            x = (x - self.x0) ./ (self.x1 - self.x0);  % (1, numPoints)
            y = (y - self.y0) ./ (self.y1 - self.y0);  % (1, numPoints)

            a = aCoef' * (cos(self.akx * x) .* cos(self.aky * y));  % (1, numPoints)
            b = bCoef' * (cos(self.bkx * x) .* sin(self.bky * y));  % (1, numPoints)
            c = cCoef' * (sin(self.ckx * x) .* cos(self.cky * y));  % (1, numPoints)
            d = dCoef' * (sin(self.dkx * x) .* sin(self.dky * y));  % (1, numPoints)
        
            fVal = reshape(a + b + c + d, originalSize);  % (*shape)
        end

        function fVal = evalBasis(self, x, y, i)
            % Evaluates a single GRF Fourier series basis function
            %
            % Parameters
            % ----------
            % x: (shape) x coordinates
            % y: (shape) y coordinates
            % i: index of basis function to evaluate
            %
            % Return
            % ------
            % fVal: (shape) basis function i evaluated (elementwise) at 
            %       (x, y)

            originalSize = size(x);
            x = reshape(x, 1, []);  % (1, numPoints)
            y = reshape(y, 1, []);  % (1, numPoints)
            
            % Transform to standard domain
            x = (x - self.x0) ./ (self.x1 - self.x0);  % (1, numPoints)
            y = (y - self.y0) ./ (self.y1 - self.y0);  % (1, numPoints)

            if i <= length(self.akx)
                fVal = cos(self.akx(i) * x) .* cos(self.aky(i) * y) ./ self.an(i);
            elseif i <= length(self.akx) + length(self.bkx)
                i = i - length(self.akx);
                fVal = cos(self.bkx(i) * x) .* sin(self.bky(i) * y) ./ self.bn(i);
            elseif i <= length(self.akx) + length(self.bkx) + length(self.ckx)
                i = i - length(self.akx) - length(self.bkx);
                fVal = sin(self.ckx(i) * x) .* cos(self.cky(i) * y) ./ self.cn(i);
            else
                i = i - length(self.akx) - length(self.bkx) - length(self.ckx);
                fVal = sin(self.dkx(i) * x) .* sin(self.dky(i) * y) ./ self.dn(i);
            end

            fVal = reshape(fVal, originalSize);  % (*shape)
        end

        function showBasis(self)
            % Shows the basis functions of the GRF
            [x, y] = meshgrid(linspace(self.x0, self.x1, 128)', linspace(self.y0, self.y1, 128)');
            
            counter = 1;
            total = self.basisDimension();
            
            % Show a basis functions
            for i = 1:length(self.akx)
                f = self.evalBasis(x, y, counter);
                surf(x, y, f);
                gx = round(self.akx(i) / (2*pi));
                gy = round(self.aky(i) / (2*pi));
                title(sprintf("a (cos(x)cos(y)) basis func %d, gx = %d, gy = %d", i, gx, gy));
                input(sprintf("Press Enter to view next (%d/%d)", counter, total));
                xlabel('x');
                ylabel('y');
                counter = counter + 1;
            end

            % Show b basis functions
            for i = 1:length(self.bkx)
                f = self.evalBasis(x, y, counter);
                surf(x, y, f);
                gx = round(self.bkx(i) / (2*pi));
                gy = round(self.bky(i) / (2*pi));
                title(sprintf("b (cos(x)sin(y)) basis func %d, gx = %d, gy = %d", i, gx, gy));
                xlabel('x');
                ylabel('y');
                input(sprintf("Press Enter to view next (%d/%d)", counter, total));
                counter = counter + 1;
            end

            % Show c basis functions
            for i = 1:length(self.ckx)
                f = self.evalBasis(x, y, counter);
                surf(x, y, f);
                gx = round(self.ckx(i) / (2*pi));
                gy = round(self.cky(i) / (2*pi));
                title(sprintf("c (sin(x)cos(y)) basis func %d, gx = %d, gy = %d", i, gx, gy));
                xlabel('x');
                ylabel('y');
                input(sprintf("Press Enter to view next (%d/%d)", counter, total));
                counter = counter + 1;
            end

            % Show d basis functions
            for i = 1:length(self.dkx)
                f = self.evalBasis(x, y, counter);
                surf(x, y, f);
                gx = round(self.dkx(i) / (2*pi));
                gy = round(self.dky(i) / (2*pi));
                title(sprintf("d (sin(x)sin(y)) basis func %d, gx = %d, gy = %d", i, gx, gy));
                xlabel('x');
                ylabel('y');
                input(sprintf("Press Enter to view next (%d/%d)", counter, total));
                counter = counter + 1;
            end
        end

        function validateOrthonormality(self)
            % Validates that the basis functions are orthonormal
            [x, y] = meshgrid(linspace(self.x0, self.x1, 129)', linspace(self.y0, self.y1, 129)');
            
            % Construct basis functions
            basis = zeros(self.basisDimension(), size(x, 1), size(x, 2));
            for i = 1:self.basisDimension()
                basis(i, 1:end, 1:end) = self.evalBasis(x, y, i);
            end

            % Evaluate inner products
            gramian = zeros(self.basisDimension(), self.basisDimension());
            for i = 1:self.basisDimension()
                for j = 1:self.basisDimension()
                    basisi = reshape(basis(i, 1:end, 1:end), [], 1);
                    basisj = reshape(basis(j, 1:end, 1:end), [], 1);
                    da = (self.x1 - self.x0) * (self.y1 - self.y0) / (128^2);
                    gramian(i, j) = (basisi' * basisj) * da;
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
            % u: function (x, y) -> u(x, y) that applies the desired
            %    function elementwise to arrays x and y
            %
            % Return
            % ------
            % uCoef: (n, 1) coefficients of the given function in the GRF
            %        basis

            [x, y] = meshgrid(linspace(self.x0, self.x1, 129)', linspace(self.x0, self.x1, 129)');
            x = reshape(x, [], 1);  % (numPoints, 1)
            y = reshape(y, [], 1);  % (numPoints, 1)
            uVal = u(x, y);  % (numPoints, 1)

            uCoef = zeros(self.basisDimension(), 1);
            for i = 1:self.basisDimension()
                basisi = self.evalBasis(x, y, i);  % (numPoints, 1)
                da = (self.x1 - self.x0) * (self.y1 - self.y0) / (128^2);
                uCoef(i) = (basisi' * uVal) * da;
            end
        end

        function n = basisDimension(self)
            % Gets the dimension of the basis used by the GRF
            %
            % Return
            % ------
            % n: Dimension of the basis

            n = length(self.akx) + length(self.bkx) + length(self.ckx) + length(self.dkx);
        end
    end
end
