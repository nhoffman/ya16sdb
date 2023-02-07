FROM python:3.11-slim-bullseye

RUN apt-get -y update && \
    apt-get -y upgrade

WORKDIR /opt/run
COPY requirements.txt .
RUN python3 -m venv /usr/local/share/venv
RUN /usr/local/share/venv/bin/python -m pip install -r requirements.txt

COPY app.py data.py entrypoint.sh filter_details.feather.gz .

EXPOSE 8000
ENTRYPOINT ["/opt/run/entrypoint.sh"]
CMD ["run"]
