FROM python:alpine

WORKDIR /tmp

RUN apk add alpine-sdk imagemagick libtool wget autoconf automake doxygen python3 python3-dev py3-pip freetype-dev
RUN ln -s /usr/bin/pip3 /usr/bin/pip
RUN ln -sf /usr/bin/pip3 /usr/local/bin/pip
RUN ln -s /usr/bin/python3 /usr/bin/python
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
RUN wget -qO- "https://www.hdfgroup.org/package/hdf5-1-10-5-tar-bz2/?wpdmdl=13570&refresh=5dde4872611c01574848626" | tar xjf -
WORKDIR /tmp/hdf5-1.10.5
RUN ./autogen.sh
RUN ./configure --enable-cxx --enable-build-mode=production --prefix=/usr
RUN make
RUN make install

RUN rm -rf /tmp/hdf5-1.10.5

## install python libraries
WORKDIR /tmp/
RUN apk add py3-numpy py-numpy-dev py3-scipy py3-jinja2
# RUN pip install --upgrade pip
RUN pip install colorama matplotlib h5py pkgconfig
