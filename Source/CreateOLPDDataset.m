
% ========================================================================
% The function CreateOLPDDataset initializes a blank .ol.h5 dataset file
% ========================================================================

% Input
% -----
% OutputFile: string, the name of the .ol.h5 file in which to create the
%             blank dataset
% GridFile: string, the name of the grid file containing the grid that will
%           be used to generate the dataset
% NumSamples: integer > 0, number of samples in the dataset 

% Output
% ------
% A .ol.h5 dataset file located at OutputFile (overwrites if file already
% exists), initialized and ready to receive NumSamples samples generated on
% the given grid.

function CreateOLPDDataset(OutputFile, x, y, uDOut, vDOut, NumSamples)
    if exist(OutputFile, "file")
        delete(OutputFile);
    end
    
    x = x';
    h5create(OutputFile, "/x/1", size(x), "Datatype", "single");
    h5write(OutputFile, "/x/1", x);
    h5writeatt(OutputFile, "/x/1", "id", int32(1));
    
    y = y';
    h5create(OutputFile, "/y/1", size(y), "Datatype", "single");
    h5write(OutputFile, "/y/1", y);
    h5writeatt(OutputFile, "/y/1", "id", int32(1));
    
    outIndices = int32(0 : (NumSamples - 1));
    h5create(OutputFile, "/u/1/indices", NumSamples, "Datatype", "int32");
    h5write(OutputFile, "/u/1/indices", outIndices);
    h5create(OutputFile, "/u/1/u", [uDOut, size(x, 2), NumSamples], "Datatype", "single");
    h5writeatt(OutputFile, "/u/1", "disc_id", int32(1));

    h5create(OutputFile, "/v/1/indices", NumSamples, "Datatype", "int32");
    h5write(OutputFile, "/v/1/indices", outIndices);
    h5create(OutputFile, "/v/1/v", [vDOut, size(y, 2), NumSamples], "Datatype", "single");
    h5writeatt(OutputFile, "/v/1", "disc_id", int32(1));
end