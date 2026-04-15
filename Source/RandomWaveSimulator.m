classdef RandomWaveSimulator < handle
    properties
       num_modes = 16
       alpha = 0.01
       beta = 1.0
       gamma = 1.9
       x
       y
       uDOut = 2
       vDOut = 2
       stream
    end

    methods
        function Setup(self, simulation, numChunks)
            self.x = simulation.xx;
            self.y = simulation.yy;
            self.stream = cell(numChunks);
            [self.stream{:}] = RandStream.create("mrg32k3a", "NumStreams", numChunks);
        end
        
        function a = Amplitude(self, i, j)
            a = self.alpha * (self.beta + i.^2 + j.^2) .^ (-self.gamma / 2);
        end

        function Generate(self, outputFile, simulation, chunkIndex, sampleIndex)
            curStream = self.stream{chunkIndex};

            % x-component of initial displacement
            voref = GRF2D(self.num_modes, @self.Amplitude, curStream);  
            vofunc = @(x, y) voref((x - Xo) / (Xn - Xo), (y - Yo) / (Yn - Yo));
            % y-component of initial displacement
            woref = GRF2D(self.num_modes, @self.Amplitude, curStream);
            wofunc = @(x, y) woref((x - Xo) / (Xn - Xo), (y - Yo) / (Yn - Yo));

            Vvofunc = @(x,y) 0.*x + 0.*y; % x-component of initial velocity
            Vwofunc = @(x,y) 0.*x + 0.*y; % y-component of initial velocity

            simulation.ImposeInitialConditions(vofunc, wofunc, Vvofunc, Vwofunc);
            u = [simulation.u, simulation.v];  % (numNodes, 2)
            
            simulation.Solve();
            v = [simulation.u, simulation.v];  % (numNodes, 2)
            
            WriteOLPDSample(outputFile, sampleIndex, u, v);
        end
    end
end