FROM python:2.7-slim

LABEL \
  author="Arnar Flatberg" \
  maintainer="arnar.flatberg@ntnu.no"

RUN \
  apt-get update && apt-get install -y --no-install-recommends \
  g++ \
  git \
  wget && \
  wget --quiet -O /opt/get-pip.py https://bootstrap.pypa.io/get-pip.py && \
  python /opt/get-pip.py && \
  rm -rf /var/lib/apt/lists/* /opt/get-pip.py

RUN pip install git+git://github.com/flatberg/MultiQC.git