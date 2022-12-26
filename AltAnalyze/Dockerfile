# syntax=docker/dockerfile:1

# pull the base image
FROM python:2

## We are building a Debian system
# PRETTY_NAME="Debian GNU/Linux 10 (buster)"
# NAME="Debian GNU/Linux"
# VERSION_ID="10"
# VERSION="10 (buster)"
# VERSION_CODENAME=buster
# ID=debian
# HOME_URL="https://www.debian.org/"
# SUPPORT_URL="https://www.debian.org/support"
# BUG_REPORT_URL="https://bugs.debian.org/"

# root is the default user
USER root  
WORKDIR /mnt

# install STAR and genome file for building index
RUN apt-get update \
    && apt-get install -y sudo xxd parallel

RUN cd /tmp \
    && wget https://github.com/alexdobin/STAR/archive/2.4.0h.tar.gz \
    && tar -xzf 2.4.0h.tar.gz \
    && sudo chown root STAR-2.4.0h \
    && cd STAR-2.4.0h/source \
    && make STAR \
    && cp STAR /usr/local/bin

RUN cd /tmp \
    && wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz \
    && gunzip gencode.v36.annotation.gtf.gz \
    && curl https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834 > GRCh38.d1.vd1.fa.tar.gz \
    && tar -xzf GRCh38.d1.vd1.fa.tar.gz 

COPY STAR_index.sh /usr/src
COPY STAR_align.sh /usr/src
COPY run.sh /usr/src

# multipath-psi
COPY requirements.txt /usr/src
COPY Hs.bed /usr/src
COPY AltAnalyzes.sh /usr/src
COPY prune.py /usr/src

RUN cd /usr/src \
    && pip install --no-cache-dir -r requirements.txt \
    && git clone https://github.com/nsalomonis/altanalyze.git \
    && cd altanalyze \
    && python AltAnalyze.py --species Hs --update Official --version EnsMart91 --additional all \
    && mv /usr/src/Hs.bed /usr/src/altanalyze/AltDatabase/EnsMart91/ensembl/Hs


ENTRYPOINT ["/usr/src/run.sh"]  # must double quote

## To build -- docker build --no-cache -t altanalyze .   (you can add -f to specify where the dockerfile is, or just build on the dir where dockerfile sit, -t image name must be lowercase)
## To push  -- docker tag altanalyze frankligy123/altanalyze:0.5.0   (first create repository on DockerHub)
## To push  -- docker push frankligy123/altanalyze:0.5.0
## To pull -- docker pull frankligy123/altanalyze
## To run  -- docker run -v $PWD:/usr/src/app/run -t frankligy123/altanalyze bam    # don't bind working dir because all the downloaded build will be wiped to be consistent with $PWD
## To run (interactive) -- docker run --rm -it --entrypoint bash frankligy123/altanalyze 
## For singularity interactive -- singularity shell --writable ./altanalyze
## To check -- docker images
## To check -- docker ps
## To check -- docker stop container_id
## To check -- docker image rm -f altanalyze
# The images help me to navigate and learn docker
# ggsashimi
# optitype
# prosit
