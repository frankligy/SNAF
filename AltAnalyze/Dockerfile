# syntax=docker/dockerfile:1

FROM python:2
WORKDIR /usr/src/app
RUN apt-get clean \
    && apt-get update \
    && apt-get install -y parallel
RUN cd /usr/src/app \
    && curl -o requirements.txt https://raw.githubusercontent.com/frankligy/SNAF/main/AltAnalyze/requirements_slim_test.txt \
    && apt-get install -y libigraph0-dev \
    && pip install python-igraph==0.7.1.post6 \
    && pip install --no-cache-dir -r requirements.txt \
    && git clone https://github.com/nsalomonis/altanalyze.git \
    && cd altanalyze \
    && python AltAnalyze.py --species Hs --update Official --version EnsMart91 --additional all \
    && cd .. \
    && curl -o Hs.bed https://raw.githubusercontent.com/frankligy/SNAF/main/AltAnalyze/Hs.bed \
    && mv ./Hs.bed ./altanalyze/AltDatabase/EnsMart91/ensembl/Hs \
    && curl -o prune.py https://raw.githubusercontent.com/frankligy/SNAF/main/AltAnalyze/prune.py \
    && chmod 777 /usr/src/app/prune.py \
    && curl -o AltAnalyze.sh https://raw.githubusercontent.com/frankligy/SNAF/main/AltAnalyze/AltAnalyze.sh \
    && chmod 777 /usr/src/app/AltAnalyze.sh \
    && chmod -R 777 /usr/src/app/altanalyze 

ENTRYPOINT ["/usr/src/app/AltAnalyze.sh"]  # must double quote






