# ---Anaconda version---
# This version works, but produces an image that is >3.1Gb
#FROM continuumio/miniconda3:latest
#RUN apt-get update && apt-get install --no-install-recommends -y \
#    build-essential \
#    gcc \
#    git \
#    && apt-get clean \
#    && rm -rf /tmp/downloaded_packages/* \
#    && rm -rf /var/lib/apt/lists/*
#
#COPY environment.yml /tmp/environment.yml
#RUN conda env create -f /tmp/environment.yml && \
#    conda clean --all -y && \
#    echo "source activate merrycrispr" > ~/.bashrc
#ENV VIRTUAL_ENV=/opt/conda/envs/merrycrispr/
#ENV PATH="$VIRTUAL_ENV/bin:$PATH"
#
#COPY . /tmp/merrycrispr
#RUN pip install git+https://github.com/lmjohns3/theanets && \
#    pip install git+https://github.com/milescsmith/Azimuth && \
#    pip install --no-cache-dir /tmp/merrycrispr && \
#    rm -fr ~/.cache/pip
#
#CMD merrycrispr
# ---Alpine version---
#FROM python:3.6.8-alpine3.10
# Not using the official python version because there is something that keeps numpy and scipy from installing
# the below image has fixed whatever the problems are
# This Alpine version is only 888 Mb
FROM frolvlad/alpine-python-machinelearning:latest

RUN apk --no-cache --update-cache add curl gcc gfortran git build-base freetype-dev libpng-dev openblas-dev python3-dev && \
    rm -rf /var/cache/apk/*
RUN ln -s /usr/include/locale.h /usr/include/xlocale.h

COPY requirements.txt ./
COPY . /tmp/merrycrispr
RUN pip install --no-cache-dir -r requirements.txt
# two separate lines here to avoid installing requirements twice
RUN pip install --no-cache-dir /tmp/merrycrispr

WORKDIR /root
RUN curl -o bowtie.zip -L https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.2/bowtie-1.2.2-linux-x86_64.zip/download && \
    unzip bowtie.zip && \
    mv bowtie-1.2.2-linux-x86_64 /usr/local/bowtie && \
    rm -rf bowtie.zip

CMD merrycrispr

# ---Debian version---
# Also works, but still ~1 Gb.
#FROM python:3.6.8-stretch
#RUN apt-get update && apt-get install --no-install-recommends -y \
#    build-essential \
#    gcc \
#    git \
#    bowtie \
#    && apt-get clean \
#    && rm -rf /tmp/downloaded_packages/* \
#    && rm -rf /var/lib/apt/lists/*
#COPY requirements.txt ./
#COPY . /tmp/merrycrispr
#RUN pip install --no-cache-dir -r requirements.txt && \
#    pip install --no-cache-dir /tmp/merrycrispr
#CMD merrycrispr