classdef RandomCrack1Simulator
    properties
        numModes = 3
        alpha = 0.1
        beta = 1.0
        gamma = 1.0
        x
        y
        dy
        uDOut = 2;
        vDOut = 1;
        stream
    end

    methods
        function Setup(self, simulator, numChunks, dy)
            self.x = simulator.xx;
            self.y = simulator.yy;
            self.dy = dy;
            self.stream = cell(numChunks);
            [self.stream{:}] = RandStream.create("mrg32k3a", "NumStreams", numChunks);
        end

        function Generate(self, outputFile, simulation, chunkIndex, sampleIndex)
            curStream = self.stream{chunkIndex};

            modtop_n = GRF1D(3, @(i) alpha ./ ((beta + i.^2).^(gamma/2)));
            modbot_n = GRF1D(3, @(i) alpha ./ ((beta + i.^2).^(gamma/2)));
            modtop = @(x) 1 + modtop_n((x - Xo) / (Xn - Xo));
            modbot = @(x) 1 + modbot_n((x - Xo) / (Xn - Xo));
            simulation.bwfunc = @(x,y,t) ( (y > Yn - self.dy) .* modtop(x) + (y < Yo + self.dy) .* modbot(x) ).*sigma.*sign(y)/self.dy; % y-component of body force density
            
            u = []
            simulation.Solve();
        end
    end
end