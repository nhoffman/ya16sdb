FROM python:3.9

RUN apt-get update && apt-get install --assume-yes --no-install-recommends \
    ca-certificates git wget

# install deenurp
RUN git clone --branch 070-python3 https://github.com/fhcrc/deenurp.git /usr/local/share/deenurp && \
    cd /usr/local/share/deenurp/ && \
    PYTHON=/usr/local/bin/python3 \
    DEENURP=/usr/local/share/deenurp/ \
    bin/bootstrap.sh /usr/local/

# install taxtastic
RUN git clone https://github.com/fhcrc/taxtastic.git /usr/local/share/taxtastic && \
    pip3 install /usr/local/share/taxtastic/

# install ya16sdb pipeline
RUN git clone https://github.com/nhoffman/ya16sdb.git /usr/local/share/ya16sdb && \
    cd /usr/local/share/ya16sdb && bin/bootstrap.sh /usr/local/

ADD SConstruct /usr/local/share/ya16sdb/

ENTRYPOINT scons --file /usr/local/share/ya16sdb/SConstruct test=yes

CMD ["--dry-run"]
