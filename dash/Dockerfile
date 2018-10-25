FROM ubuntu:18.04
MAINTAINER crosenth@gmail.com
RUN apt-get update && apt-get install --assume-yes python3 python3-pip
COPY ./ /dash/
RUN pip3 install --requirement /dash/requirements.txt
CMD ["gunicorn", "app:app.server", "--config", "/dash/config.py"]
EXPOSE 443
