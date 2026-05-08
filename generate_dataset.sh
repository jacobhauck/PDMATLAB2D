# Configuration variables
OUTPUT_FILE=$1
GRID_FILE=GridFile1.mat
DATASET_SIZE=$2
NUM_CHUNKS=1
RANDOM_SEED=$4
CHUNK_SIZE=$(( (DATASET_SIZE + NUM_CHUNKS - 1) / NUM_CHUNKS ))
INPUT_DECK=$3

echo "Creating dataset file with chunks of size " $CHUNK_SIZE
matlab -batch "CreateOLPDDataset('$OUTPUT_FILE', '$GRID_FILE', $DATASET_SIZE)"

# Generate dataset in parallel chunks
i=1
total_data=0
while (( total_data < DATASET_SIZE )); do
    echo "Dispatching chunk " $i
    matlab -batch "OLPDMATLAB2D('$INPUT_DECK', '$OUTPUT_FILE', $CHUNK_SIZE, $i, $RANDOM_SEED)" &
    i=$((i+1))
    total_data=$((total_data + CHUNK_SIZE))
done
wait

echo "Merging dataset chunks"
matlab -batch "MergeOLPDChunks('$OUTPUT_FILE', $CHUNK_SIZE)"

echo "Removing temporary chunk files"
rm $OUTPUT_FILE.*.ol.h5

echo "Finished"

