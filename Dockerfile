FROM python:alpine

WORKDIR /tmp

RUN apk add alpine-sdk imagemagick libtool wget autoconf automake doxygen python3 python3-dev py3-pip freetype-dev linux-headers build-base libexecinfo-dev
RUN ln -sf /usr/bin/pip3 /usr/bin/pip
RUN ln -sf /usr/bin/pip3 /usr/local/bin/pip
RUN ln -sf /usr/bin/python3 /usr/bin/python
RUN ln -sf /usr/bin/python3 /usr/local/bin/python
RUN ln -sf /usr/bin/python3 /usr/local/bin/python3

## install pandoc
WORKDIR /tmp/
RUN wget -qO- https://github.com/jgm/pandoc/releases/download/2.8/pandoc-2.8-linux-amd64.tar.gz | tar xzf -
WORKDIR /tmp/pandoc-2.8
RUN for i in `find . -type f`; do D=`dirname $i`; echo mv $i /usr/${D}; done
# WORKDIR /tmp/
RUN rm -rf /tmp/pandoc-2.8

## build and install HDF5 library
WORKDIR /tmp/
RUN wget -qO- "https://www.hdfgroup.org/package/hdf5-1-12-0-tar-bz2/?wpdmdl=14584&refresh=601a6b608b5841612344160" | tar xjf -
WORKDIR /tmp/hdf5-1.12.0
RUN ./autogen.sh
RUN ./configure --enable-build-mode=production --prefix=/usr
RUN make
RUN make install

RUN rm -rf /tmp/hdf5-1.12.0

## install python libraries
WORKDIR /tmp/
RUN apk add py3-wheel py3-numpy py3-numpy-dev py3-scipy py3-jinja2 py3-pygments py3-colorama py3-matplotlib
RUN pip install --upgrade pip
COPY lib/python/picopic/setup.py setup_lib.py
COPY test/functional/setup.py setup_test.py
RUN python setup_lib.py install
RUN python setup_test.py install
RUN rm -f setup_lib.py setup_test.py

# RUN pip install colorama h5py pkgconfig

## required for stacktracing in loguru
RUN echo "export LIBRARIES=execinfo" >> /etc/profile