FROM continuumio/miniconda3:latest
RUN apt-get update && apt-get install -y git
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    git \
    wget \
    && apt-get clean \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*

ADD environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && \
    conda clean --all -y
RUN echo "source activate $(head -1 /tmp/environment.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /tmp/environment.yml | cut -d' ' -f2)/bin:$PATH

ADD . /tmp/merrycrispr
RUN pip install --no-cache-dir /tmp/merrycrispr && \
    rm -fr ~/.cache/pip

ENTRYPOINT [ "/opt/conda/envs/merrycrispr/bin", "--" ]
CMD [ "python merrycrispr" ]