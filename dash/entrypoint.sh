#!/bin/bash
set -e

source /usr/local/share/venv/bin/activate

if [[ "$1" == 'run' ]]; then
    exec gunicorn app:server \
         -b 0.0.0.0:8000 \
         -t 240 \
         --capture-output --log-file -
else
    exec "$@"
fi
