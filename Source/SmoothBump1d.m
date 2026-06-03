classdef SmoothBump1d
    properties
        n
        width
        height
        coef = []
    end

    methods
        function self = SmoothBump1d(n, width, height)
            self.n = n;
            self.width = width;
            self.height = height;
            if n > 0
                matrix = zeros(n, n);
                nVals = n : (2*n - 1);
                for row = 1:n
                    matrix(row, :) = falling(nVals, row - 1);
                end
                a = self.width / 2;
                rhs = (height / (a^n)) .* ([1; zeros(n-1, 1)]);
                c_scaled = matrix \ rhs;
                self.coef = c_scaled .* (a .^ (-(0 : (n-1))))';
            end
        end
        
        function result = evalLeft(self, x)
            cur_power = ones(size(x));
            result = zeros(size(x));
            d = x + self.width/2;
            for j = 1:self.n
                result = result + self.coef(j) * cur_power;
                cur_power = cur_power .* d;
            end
            result = result .* cur_power;
        end

        function result = evaluate(self, x)
            if self.n > 0
                resultLeft = (x < 0) .* (x >= -self.width/2) .* self.evalLeft(x);
                resultRight = (x >= 0) .* (x <= self.width/2) .* self.evalLeft(-x);
                result = resultLeft + resultRight;
            else
                result = self.height * (abs(x) < self.width/4);
            end
        end
    end
end

function result = falling(x, p)
    result = ones(size(x));
    for j = 0 : (p-1)
        result = result .* (x - j);
    end
end