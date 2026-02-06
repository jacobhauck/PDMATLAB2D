
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

function OLPDMATLAB2D(InputDeck, DatasetFile, ChunkSize, ChunkIndex)
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

        cd('InputFiles/')
        run(InputDeck);
        cd('..');

        chunkFile = sprintf("../%s.%d.ol.h5", DatasetFile, ChunkIndex);
        fprintf("Creating temporary chunk dataset at %s\n", chunkFile);
        cd('Source');
        CreateOLPDDataset(chunkFile, GridFile, effectiveChunkSize);
        cd('..');
        chunkFile = sprintf("%s.%d.ol.h5", DatasetFile, ChunkIndex);
    end
    
    for index = 1:effectiveChunkSize
        % Close all figures
        close all
    
        % Clear variables from memory
        clearvars a* -except InputDeck
    
        % Clear command window
        clc
    
        % Print simulation title
        fprintf(' ================================================================ \n')
        fprintf('                          %s \n',InputDeck)
        fprintf(' ================================================================ \n\n')
    
        % Read input deck for specific simulation problem
        cd('InputFiles/')
        run(InputDeck);
    
        % Run Main script
        cd('../Source/')
        Main
     
        cd .. 
        
        uOut = [vInitial, wInitial];  % (numNodes, 2)
        vOut = [v, w];  % (numNodes, 2)
        WriteOLPDSample(chunkFile, index, uOut, vOut);
    end
end