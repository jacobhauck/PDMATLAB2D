
% ========================================================================
% Copyright (c) 2022 by Oak Ridge National Laboratory                      
% All rights reserved.                                                     
%                                                                           
% This file is part of PDMATLAB2D. PDMATLAB2D is distributed under a           
% BSD 3-clause license. For the licensing terms see the LICENSE file in    
% the top-level directory.                                                 
%                                                                          
% SPDX-License-Identifier: BSD-3-Clause                                    
% ========================================================================

function OLPDMATLAB2D(InputDeck, DatasetFile, ChunkSize, ChunkIndex, Seed)
    numSamples = length(h5read(DatasetFile, "/u/1/indices"));
    
    effectiveChunkSize = numSamples;
    if nargin == 2
        chunkFile = DatasetFile;
    else
        offset = ChunkSize * (ChunkIndex - 1);
        if offset >= numSamples
            error("Chunk index too large for dataset size");
        end
        effectiveChunkSize = min(ChunkSize, numSamples - offset);

        if effectiveChunkSize < numSamples
            fprintf("Chunk %d: initializing independent stream of random numbers\n", ChunkIndex);
            numChunks = idivide(int64(numSamples + ChunkSize - 1), int64(ChunkSize));
            streams = cell(1, numChunks);
            % Create full set of streams; with a fixed seed, this is
            % deterministic across all chunks. Then take the stream for the
            % current ChunkIndex
            [streams{:}] = RandStream.create("mrg32k3a", "NumStreams", numChunks, "Seed", Seed);
            RandStream.setGlobalStream(streams{ChunkIndex});
        end

        run(['InputFiles/', InputDeck]);
        fprintf("Chunk %d: using grid file %s\n", ChunkIndex, GridFile);
        chunkFile = sprintf("%s.%d.ol.h5", DatasetFile, ChunkIndex);
        fprintf("Chunk %d: creating temporary chunk dataset at %s\n", ChunkIndex, chunkFile);
        CreateOLPDDataset(chunkFile, GridFile, effectiveChunkSize);
    end
    
    for index = 1:effectiveChunkSize
        fprintf("Chunk %d: running sample %d/%d\n", ChunkIndex, index, effectiveChunkSize);
        
        % Generate initial condition
        run(['InputFiles/', InputDeck]);

        % Run Main script
        run('Source/Main');
        
        % Save initial and final states
        uOut = [vInitial, wInitial];  % (numNodes, 2)
        vOut = [v, w];  % (numNodes, 2)
        WriteOLPDSample(chunkFile, index, uOut, vOut);
    end
end