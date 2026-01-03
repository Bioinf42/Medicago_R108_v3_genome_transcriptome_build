#!/usr/bin/env bash
set -euo pipefail

IMAGE="${IMAGE:-txpipe:0.3}"
ENV_NAME="${ENV_NAME:-transcriptomeDockerPipe}"

mkdir -p raw_data reference results logs results/fastq_untrimmed results/trimmed

docker run --rm -it \
  --user "$(id -u)":"$(id -g)" \
  -e HOME=/tmp \
  -e XDG_CACHE_HOME=/tmp/.cache \
  -e MAMBA_USER=true \
  -e STEP \
  -e THREADS \
  -v "$(pwd):/work" \
  --entrypoint micromamba \
  "$IMAGE" \
  run -n "$ENV_NAME" /work/transcriptome_build_pipe.sh

