#!/bin/bash

# ./jupiter_launch.sh

# [INFO] Setting up environment variables
CONTAINER_NAME="jupiter_images_Deqtectseq"
IMAGE_NAME="detectseqpipe:6"
HOST_PATH="/30tb/home/nuvobea/pmat_and_mpmat/CH12F3/Scripts_DepIndep/jupiter_images"
RESULTS_PATH="/30tb/home/nuvobea/pmat_and_mpmat/CH12F3/DepIndep_dataset"
HOST_PORT=8001

# [INFO] Checking if container already exists
if [ "$(docker ps -aq -f name=${CONTAINER_NAME})" ]; then
    echo "[WARNING] Removing existing container: ${CONTAINER_NAME}"
    docker rm -f ${CONTAINER_NAME}
fi

# [INFO] Starting the Docker container
# Note: Removed python3 -m webbrowser as it won't work over SSH/VS Code
docker run -d \
    --name "${CONTAINER_NAME}" \
    --platform linux/amd64 \
    -v "${HOST_PATH}":/sharedFolder \
    -v "${RESULTS_PATH}":/sharedFolder/resultsFolder \
    --memory-reservation=10g \
    --privileged=true \
    -p "${HOST_PORT}":8888 \
    "${IMAGE_NAME}"

if [ $? -eq 0 ]; then
    echo "[INFO] Container started successfully."
    echo "[INFO] Access Jupyter Lab locally at: http://localhost:${HOST_PORT}"
else
    echo "[ERROR] Failed to start the container."
    exit 1
fi

# [INFO] To access Jupyter Lab remotely, set up SSH tunneling with the following command:
# ssh -L 8001:localhost:8001 utente@indirizzo-del-server
# ssh -L 8001:localhost:8001 nuvoba@130.192.212.153