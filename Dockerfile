FROM python:3.13-bookworm
ARG YA16SDB_VERSION
ENV \
PIP_ROOT_USER_ACTION=ignore \
SCONSFLAGS="--file /usr/local/share/ya16sdb/SConstruct" \
YA16SDB_VERSION=${YA16SDB_VERSION}
RUN apt-get update && \
apt-get upgrade -y && \
apt-get remove -y curl && \
apt-get install -y ca-certificates wget
WORKDIR /usr/local/share/ya16sdb/
COPY requirements.txt SConstruct ncbi.conf ./
COPY testfiles/ ./testfiles/
COPY tests/ ./tests/
COPY bin/ ./bin/
COPY data/ ./data/
RUN bin/bootstrap.sh /usr/local/
RUN find . -type f -exec chmod 644 {} \; && \
find . -type d -exec chmod 755 {} \; && \
find ./bin/ -type f -exec chmod 755 {} \;
CMD ["scons", "--dry-run"]
