# Custom Modkit container with Python dependencies for transcriptome acceleration
FROM ontresearch/modkit:shaa7bf2b62946eeb7646b9b9d60b892edfc3b3a52c

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install micromamba
RUN curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

# Copy conda environment file
COPY conda.yml /opt/conda.yml

# Create environment from YAML file, clean cache, and generate lock file
RUN micromamba install -y -n base -f /opt/conda.yml \
    && micromamba install -y -n base "setuptools>=78.1.1" \
    && micromamba clean -afy \
    && micromamba env export -n base --explicit > /opt/conda.lock

# Activate the environment by default
SHELL ["/usr/local/bin/micromamba", "run", "-n", "base", "/bin/bash", "-c"]

# Initialize micromamba for interactive shells
RUN micromamba shell init --shell bash \
    && echo "micromamba activate base" >> ~/.bashrc

# Set the working directory used by Nextflow processes
WORKDIR /workspace

# Set default command to run bash within the activated environment
CMD ["micromamba", "run", "-n", "base", "bash"]

