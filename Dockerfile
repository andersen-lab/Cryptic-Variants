FROM --platform=linux/x86_64 mambaorg/micromamba:0.25.1

LABEL image.author.name "Dylan Pilz"
LABEL image.author.email "dpilz@scripps.edu"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

RUN micromamba create -n cryptic-variants

RUN micromamba install -y -n cryptic-variants -f /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/cryptic-variants/bin:$PATH