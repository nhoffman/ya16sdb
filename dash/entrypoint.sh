#!/bin/bash
set -e

source /usr/local/share/venv/bin/activate

if [[ "$1" == 'run' ]]; then
    exec gunicorn app:server \
         -b 0.0.0.0:8000 \
         --timeout 600 \
         --capture-output --log-file - & (dmesg | grep gunicorn)
else
    exec "$@"
fi
