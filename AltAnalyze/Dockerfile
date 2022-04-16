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
    && curl -o AltAnalyze.sh https://raw.githubusercontent.com/frankligy/SNAF/main/AltAnalyze/AltAnalyze.sh
ENTRYPOINT ["AltAnalyze.sh"]

## To build -- docker build -t AltAnalyze .

## To pull -- docker pull frankligy123/AltAnalyze

## To run  -- docker run -v $PWD:/usr/src/app -t frankligy123/AltAnalyze ./bam
