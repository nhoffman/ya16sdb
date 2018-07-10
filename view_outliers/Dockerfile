FROM ubuntu:18.04
MAINTAINER crosenth@gmail.com
RUN apt-get update && apt-get install --assume-yes python3 python3-pip
COPY ./ /view_outliers/
RUN pip3 install --requirement /view_outliers/requirements.txt
CMD ["bokeh", "serve", "--allow-websocket-origin", "localhost:8006", "--port", "8005", "/view_outliers"]
EXPOSE 8005
