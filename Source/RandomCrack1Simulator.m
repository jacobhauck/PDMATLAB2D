classdef RandomCrack1Simulator
    properties
        x
        y
        uDOut = 2;
        vDOut = 1;
        stream
    end

    methods
        function Setup(self, simulator, numChunks)
            self.x = simulator.xx;
            self.y = simulator.yy;
            self.stream = cell(numChunks);
            [self.stream{:}] = RandStream.create("mrg32k3a", "NumStreams", numChunks);
        end

        function Generate(self, outputFile, simulation, chunkIndex, sampleIndex)
            curStream = self.stream{chunkIndex};
        end
    end
end