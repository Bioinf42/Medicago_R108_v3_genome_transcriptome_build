FROM mambaorg/micromamba:latest

# Copy your environment recipe into the image
COPY docker/env.yaml /tmp/env.yaml

# Create the environment described in env.yaml
RUN micromamba create -y -n transcriptomeDockerPipe -f /tmp/env.yaml && micromamba clean -a -y

# Make the tools in that environment available by default
ENV PATH=/opt/conda/envs/transcriptomeDockerPipe/bin:$PATH

# Pick a default folder inside the container
WORKDIR /work

# Default command: show FastQC version (a simple test)
ENTRYPOINT ["micromamba", "run", "-n", "transcriptomeDockerPipe"]
CMD ["fastqc", "--version"]

