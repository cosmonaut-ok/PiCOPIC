#!/bin/sh

DISP=42
if [ ! -z ${1} ]; then
  DISP=${1}
fi

Xvfb :${DISP} -screen 0 1500x990x16 &

echo "export DISPLAY=:${DISP}"
