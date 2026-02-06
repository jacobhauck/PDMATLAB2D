
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
    
    if nargin == 2
        ChunkSize = numSamples;
        ChunkIndex = 1;
    end

    offset = ChunkSize * (ChunkIndex - 1);
    if offset >= numSamples
        error("Chunk index too large for dataset size");
    end

    effectiveChunkSize = min(ChunkSize, numSamples - offset);
    
    for irel = 1:effectiveChunkSize
        index = offset + irel;

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
        WriteOLPDSample(DatasetFile, index, uOut, vOut);
    end
end