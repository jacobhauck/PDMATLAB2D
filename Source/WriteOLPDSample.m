
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
    uNumNodes = size(u, 1);
    uDOut = size(u, 2);
    uData = single(reshape(u', uDOut, numNodes, 1));  % (uDOut, uNumNodes, 1)
    h5write(DatasetFile, "/u/1/u", uData, [1, 1, Index], [uDOut, uNumNodes, 1]);
    
    vNumNodes = size(u, 1);
    vDOut = size(v, 2);
    vData = single(reshape(v', vDOut, numNodes, 1));  % (vDOut, vNumNodes, 1)
    h5write(DatasetFile, "/v/1/v", vData, [1, 1, Index], [vDOut, vNumNodes, 1]);
end