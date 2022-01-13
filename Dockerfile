FROM python:3.9

RUN apt-get update && apt-get install --assume-yes --no-install-recommends ca-certificates git wget

# Add files
RUN git clone https://github.com/nhoffman/ya16sdb.git /usr/local/share/ya16sdb && \
    cd /usr/local/share/ya16sdb && \
    bin/bootstrap.sh /usr/local/

# ADD SConstruct /usr/local/share/ya16sdb/

RUN git clone https://github.com/fhcrc/taxtastic.git /usr/local/share/taxtastic && \
    pip3 install /usr/local/share/taxtastic/

CMD ["scons", "--dry-run", "--file", "/usr/local/share/ya16sdb/SConstruct"]
