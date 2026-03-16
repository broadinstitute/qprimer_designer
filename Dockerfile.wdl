FROM python:3.12-slim-bookworm

LABEL org.opencontainers.image.source="https://github.com/broadinstitute/qprimer_designer"
LABEL org.opencontainers.image.description="WDL runner for the qPrimer Designer pipeline"
LABEL org.opencontainers.image.licenses="MIT"

# Install miniwdl and its dependencies
RUN pip install --no-cache-dir miniwdl>=1.11

# Copy WDL workflow files
COPY workflows/*.wdl /opt/qprimer/workflows/
COPY workflows/*.json /opt/qprimer/workflows/
COPY workflows/params.txt.template /opt/qprimer/workflows/

# Default miniwdl config: use Docker backend, pull task image automatically
RUN mkdir -p /root/.config && \
    printf '[scheduler]\ncall_concurrency = 16\n' > /root/.config/miniwdl.cfg

# Verify miniwdl installation and WDL syntax
RUN miniwdl check /opt/qprimer/workflows/qprimer.wdl

ENV QPRIMER_WDL_DIR=/opt/qprimer/workflows

WORKDIR /data

ENTRYPOINT ["miniwdl", "run", "/opt/qprimer/workflows/qprimer.wdl"]
CMD ["--help"]
