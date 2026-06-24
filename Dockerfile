FROM mambaorg/micromamba:2.5.0-debian12-slim

LABEL org.opencontainers.image.source="https://github.com/broadinstitute/qprimer_designer"
LABEL org.opencontainers.image.description="ML-guided qPCR primer design with off-target minimization"
LABEL org.opencontainers.image.licenses="MIT"

# Layer 1: Conda environment (cached unless environment.yml or TORCH_VARIANT changes)
#
# TORCH_VARIANT selects the PyTorch build:
#   gpu (default) -> CUDA-enabled build; used for the multi-arch GHCR image that
#                    powers training and the Terra/batch CLI (which can use GPUs).
#   cpu           -> CPU-only build (adds the `cpuonly` package); used for the slim
#                    amd64 GAR image that runs the Streamlit web app on Cloud Run
#                    (no GPU there). `cpuonly` is applied in the SAME solve so the
#                    CUDA libraries are never written to a layer -- removing them in
#                    a later layer would not shrink the image.
ARG TORCH_VARIANT=gpu
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN EXTRA=""; [ "$TORCH_VARIANT" = "cpu" ] && EXTRA="cpuonly"; \
    micromamba install -y -n base -f /tmp/environment.yml $EXTRA && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH="/opt/conda/bin:$PATH"

# Layer 2: Install package (cached unless pyproject.toml, README.md, or src/ changes)
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md /app/
COPY --chown=$MAMBA_USER:$MAMBA_USER src/ /app/src/
WORKDIR /app
RUN pip install --no-cache-dir ".[gui]"

# Layer 3: Copy remaining files (workflows, training, CLAUDE.md, LICENSE, etc.)
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

# Verify installation
RUN qprimer --help && \
    qprimer generate --help && \
    RNAduplex --version && \
    bowtie2 --version

WORKDIR /data

# Streamlit configuration for Cloud Run
ENV STREAMLIT_SERVER_PORT=8080
ENV STREAMLIT_SERVER_ADDRESS=0.0.0.0
ENV STREAMLIT_SERVER_HEADLESS=true
ENV STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

EXPOSE 8080

CMD ["sh", "-c", "streamlit run --server.port=${PORT:-8080} --server.address=0.0.0.0 /app/gui/app.py"]
