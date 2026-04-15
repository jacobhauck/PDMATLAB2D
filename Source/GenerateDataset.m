function GenerateDataset(outputFile, datasetSize, numChunks, baseSimulation, simulator)
    minChunkSize = floor(datasetSize / numChunks);
    numExtras = mod(datasetSize, numChunks);

    fprintf("Creating dataset file at: %s\n", outputFile);
    CreateOLPDDataset(outputFile, simulator.x, simulator.y, simulator.uDOut, simulator.vDOut, datasetSize);
    
    fprintf("Generating data with %d chunks of size %d\n", numChunks, minChunkSize);
    parfor i=1:numChunks
        ChunkSize = minChunkSize + (i <= numExtras);
        chunkFile = sprintf("%s-%d.ol.h5", outputFile, i);
        CreateOLPDDataset( ...
            chunkFile, ...
            simulator.x, simulator.y, simulator.uDOut, simulator.vDOut, ...
            datasetSize ...
        ); %#ok<PFBNS>
        
        for j=1:ChunkSize
            simulation = copy(baseSimulation);
            simulator.Generate(chunkFile, simulation, i, j);
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