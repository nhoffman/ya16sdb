FROM python:3.9

RUN apt-get update && apt-get upgrade -y && apt-get install --assume-yes --no-install-recommends \
    ca-certificates git wget

ADD requirements.txt /usr/local/share/ya16sdb/
ADD bin/bootstrap.sh /usr/local/share/ya16sdb/bin/

RUN cd /usr/local/share/ya16sdb && bin/bootstrap.sh /usr/local/

ADD .git/ /usr/local/share/ya16sdb/.git/
ADD data/ /usr/local/share/ya16sdb/data/
ADD bin/ /usr/local/share/ya16sdb/bin/
ADD SConstruct ncbi.conf /usr/local/share/ya16sdb/

ENV SCONSFLAGS="--file /usr/local/share/ya16sdb/SConstruct"
