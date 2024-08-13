FROM python:3.11-bookworm

ENV PIP_ROOT_USER_ACTION=ignore

RUN apt-get update && apt-get upgrade -y && \
apt-get install --assume-yes ca-certificates git wget

COPY requirements.txt /usr/local/share/ya16sdb/
COPY bin/bootstrap.sh /usr/local/share/ya16sdb/bin/

WORKDIR /usr/local/share/ya16sdb/
RUN ["/bin/bash", "-c", "bin/bootstrap.sh /usr/local/"]

COPY .git/ /usr/local/share/ya16sdb/.git/
COPY data/ /usr/local/share/ya16sdb/data/
COPY bin/ /usr/local/share/ya16sdb/bin/
COPY SConstruct ncbi.conf /usr/local/share/ya16sdb/

RUN find /usr/local/share/ya16sdb/ -type f -exec chmod 644 {} \; && \
find /usr/local/share/ya16sdb/ -type d -exec chmod 755 {} \; && \
find /usr/local/share/ya16sdb/bin/ -type f -exec chmod 755 {} \;

ENV SCONSFLAGS="--file /usr/local/share/ya16sdb/SConstruct"

CMD ["scons", "--dry-run"]
