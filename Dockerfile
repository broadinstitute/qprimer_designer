FROM mambaorg/micromamba:2.0.0-ubuntu24.04

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH="/opt/conda/bin:$PATH"

COPY --chown=$MAMBA_USER:$MAMBA_USER . /app
WORKDIR /app

RUN pip install --no-cache-dir .

# Verify installation
RUN qprimer --help && \
    qprimer generate --help && \
    RNAduplex --version && \
    bowtie2 --version | head -1

WORKDIR /data
