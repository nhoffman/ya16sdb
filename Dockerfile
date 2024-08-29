FROM python:3.11-bookworm
ARG MEFETCH_API_KEY MEFETCH_EMAIL YA16SDB_VERSION
ENV \
MEFETCH_API_KEY=${MEFETCH_API_KEY} \
MEFETCH_EMAIL=${MEFETCH_EMAIL} \
PIP_ROOT_USER_ACTION=ignore \
SCONSFLAGS="--file /usr/local/share/ya16sdb/SConstruct" \
YA16SDB_VERSION=${YA16SDB_VERSION}
RUN apt-get update && apt-get upgrade -y && apt-get install -y ca-certificates wget
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
