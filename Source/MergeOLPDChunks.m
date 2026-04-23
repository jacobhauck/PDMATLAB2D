function MergeOLPDChunks(DatasetFile, NumChunks)
    start = [1, 1, 1];
    for chunkIndex = 1:NumChunks
        chunkFile = sprintf("%s-%d.ol.h5", DatasetFile, chunkIndex);
        if ~exist(chunkFile, "file")
            error(sprintf("Cannot merge incomplete set of chunks, missing: %s", chunkFile));
        end

        u = h5read(chunkFile, "/u/1/u");
        v = h5read(chunkFile, "/v/1/v");
        
        if length(size(u)) == 2
            sizeU = [size(u), 1];
        else
            sizeU = size(u);
        end

        if length(size(v)) == 2
            sizeV = [size(v), 1];
        else
            sizeV = size(v);
        end
        h5write(DatasetFile, "/u/1/u", u, start, sizeU);
        h5write(DatasetFile, "/v/1/v", v, start, sizeV);
        
        start(3) = start(3) + sizeU(3);
    end
end
