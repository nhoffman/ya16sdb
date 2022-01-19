FROM python:3.9

RUN apt-get update && apt-get install --assume-yes --no-install-recommends \
    ca-certificates git wget

# install ya16sdb pipeline
RUN git clone https://github.com/nhoffman/ya16sdb.git /usr/local/share/ya16sdb && \
    cd /usr/local/share/ya16sdb && bin/bootstrap.sh /usr/local/

ENTRYPOINT scons --file /usr/local/share/ya16sdb/SConstruct

CMD ["--dry-run"]
