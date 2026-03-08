FROM mambaorg/micromamba:1.5.10

ARG NEXTFLOW_VERSION=25.10.4

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    bzip2 \
    ca-certificates \
    coreutils \
    curl \
    gzip \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY --chmod=0755 entrypoint.sh /usr/local/bin/rnaseq-pipeline

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/rnaseq.yml /tmp/envs/rnaseq.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/cufflinks.yml /tmp/envs/cufflinks.yml

USER $MAMBA_USER

RUN micromamba config set remote_connect_timeout_secs 30 && \
    micromamba config set remote_max_retries 10 && \
    micromamba config set remote_backoff_factor 2 && \
    micromamba create --retry-clean-cache -y -n rnaseq -f /tmp/envs/rnaseq.yml && \
    micromamba create --retry-clean-cache -y -n cufflinks -f /tmp/envs/cufflinks.yml && \
    micromamba clean --all --yes

ENV PATH=/opt/conda/envs/rnaseq/bin:/usr/local/bin:${PATH}

RUN NXF_VER=${NEXTFLOW_VERSION} bash -c 'curl -fsSL https://get.nextflow.io | bash' && \
    mv nextflow /opt/conda/envs/rnaseq/bin/nextflow && \
    chmod 0755 /opt/conda/envs/rnaseq/bin/nextflow

WORKDIR /workspace

ENTRYPOINT ["/usr/local/bin/rnaseq-pipeline"]
CMD ["help"]
