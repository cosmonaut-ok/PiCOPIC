#!/bin/sh

DISP=42
if [ ! -z ${1} ]; then
  DISP=${1}
fi

Xnest :${DISP} -ac -geometry 1500x990 &

echo "export DISPLAY=:${DISP}"
