
% ========================================================================
% The function WriteOLPDSample writes a single OL dataset sample to a
% .ol.h5 dataset file
% ========================================================================

% Input
% -----
% DatasetFile: string, the name of the .ol.h5 file in which to write the
%              sample
% Index: integer, index of the sample in the dataset
% u: array of shape (numNodes, 2), the initial displacement field
% v: array of shape (numNodes, 2), the final displacement field

% Output
% ------
% Writes the input--output pair (u, v) to the given index in the given
% .ol.h5 dataset file (from 1 to total number of samples; not the internal 
% index, which is separately handled by CreateOLPDDataset)

function WriteOLPDSample(DatasetFile, Index, u, v)
    numNodes = size(u, 1);
    if numNodes ~= size(v, 1)
        error("Invalid data: numNodes is not the same for u and v");
    end

    uData = single(reshape(u', 2, numNodes, 1));  % (2, numNodes, 1)
    h5write(DatasetFile, "/u/1/u", uData, [1, 1, Index], [2, numNodes, 1]);

    vData = single(reshape(v', 2, numNodes, 1));  % (2, numNodes, 1)
    h5write(DatasetFile, "/v/1/v", vData, [1, 1, Index], [2, numNodes, 1]);
end