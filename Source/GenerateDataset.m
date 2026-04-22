function GenerateDataset(outputFile, datasetSize, numChunks, baseSimulation, x, y, uDOut, vDOut, generator)
    minChunkSize = floor(datasetSize / numChunks);
    numExtras = mod(datasetSize, numChunks);

    fprintf("Creating dataset file at: %s\n", outputFile);
    CreateOLPDDataset(outputFile, simulator.x, simulator.y, simulator.uDOut, simulator.vDOut, datasetSize);
    
    fprintf("Generating data with %d chunks of size %d\n", numChunks, minChunkSize);
    stream = cell(numChunks);
    [stream{:}] = RandStream.create("mrg32k3a", "NumStreams", numChunks);
    parfor i=1:numChunks
        ChunkSize = minChunkSize + (i <= numExtras);
        chunkFile = sprintf("%s-%d.ol.h5", outputFile, i);
        CreateOLPDDataset( ...
            chunkFile, ...
            x, y, uDOut, vDOut, ...
            datasetSize ...
        );
        
        for j=1:ChunkSize
            simulation = copy(baseSimulation);
            generator(simulation, chunkFile, i, j, stream{i}); %#ok<PFBNS>
        end
    end
    fprintf("Finished generating data\n");

    fprintf("Merging %d chunks into dataset file %s\n", numChunks, outputFile);
    MergeOLPDChunks(outputFile, numChunks);
    
    fprintf("Removing %d temporary files\n", numChunks);
    for i=1:numChunks
        chunkFile = sprintf("%s-%d.ol.h5", outputFile, i);
        delete(chunkFile);
    end
end