# syntax=docker/dockerfile:1

FROM python:2
WORKDIR /usr/src/app
RUN curl -o requirements.txt https://raw.githubusercontent.com/frankligy/SNAF/main/AltAnalyze/requirements_slim.txt \
    && pip install --no-cache-dir -r requirements.txt \
    && git clone https://github.com/nsalomonis/altanalyze.git \
    && cd altanalyze \
    && python AltAnalyze.py --species Hs --update Official --version EnsMart91 --additional all \
    && cd .. \
    && curl -o Hs.bed https://raw.githubusercontent.com/frankligy/SNAF/main/AltAnalyze/Hs.bed \
    && mv ./Hs.bed ./altanalyze/AltDatabase/EnsMart91/ensembl/Hs \
    && curl -o AltAnalyze.sh https://raw.githubusercontent.com/frankligy/SNAF/main/AltAnalyze/AltAnalyze.sh \
    && curl -o prune.py https://raw.githubusercontent.com/frankligy/SNAF/main/AltAnalyze/prune.py \
    && chmod 777 /usr/src/app/AltAnalyze.sh \
    && mkdir run \
    && chmod -R 777 /usr/src/app/altanalyze   
ENTRYPOINT ["/usr/src/app/AltAnalyze.sh"]  # must double quote

## To build -- docker build --no-cache -t altanalyze .

## To push  -- docker tag altanalyze frankligy123/altanalyze:0.5.0   (first create repository on DockerHub)
## To push  -- docker push frankligy123/altanalyze:0.5.0

## To pull -- docker pull frankligy123/altanalyze

## To run  -- docker run -v $PWD:/usr/src/app/run -t frankligy123/altanalyze bam    # don't bind working dir because all the downloaded build will be wiped to be consistent with $PWD

## To run (interactive) -- docker run --rm -it --entrypoint bash frankligy123/altanalyze 

## To check -- docker images
## To check -- docker ps
## To check -- docker image rm -f altanalyze
