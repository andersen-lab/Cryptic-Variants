FROM mambaorg/micromamba:0.25.1

LABEL image.author.name "Dylan Pilz"
LABEL image.author.email "dpilz@scripps.edu"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

RUN micromamba create -n freyja

RUN micromamba install -y -n freyja -f /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/freyja/bin:$PATH