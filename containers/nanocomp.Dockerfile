
FROM dhi.io/debian-base:bookworm

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    bzip2 \
    && rm -rf /var/lib/apt/lists/*

# Install micromamba
RUN mkdir -p /usr/local/bin && \
    curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

# Install nanocomp directly into base environment
RUN micromamba install -y -n base -c bioconda -c conda-forge nanocomp \
    && micromamba clean -afy

ENV PATH="/root/.local/share/mamba/bin:${PATH}"

# Activate the environment by default
SHELL ["/usr/local/bin/micromamba", "run", "-n", "base", "/bin/bash", "-c"]

# Initialize micromamba for interactive shells
RUN micromamba shell init --shell bash \
    && echo "micromamba activate base" >> ~/.bashrc

# Set default command
CMD ["micromamba", "run", "-n", "base", "bash"]

