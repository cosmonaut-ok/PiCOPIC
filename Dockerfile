FROM debian:buster-slim

WORKDIR /tmp

RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

RUN echo 'deb-src http://deb.debian.org/debian stretch main' >> /etc/apt/sources.list
RUN apt-get update -y -qq

RUN mkdir -p /usr/share/man/man1
RUN apt-get -y -qq install apt-utils build-essential git imagemagick pandoc pandoc-citeproc libtool libtool-bin wget autoconf doxygen pandoc-citeproc python3-numpy python3-colorama python3-jinja2 python3-scipy python3-matplotlib || true

# because mainstream openjdk does not want to be installed w/o some (pseudo)graphics
# RUN apt-get -y -qq install default-jdk-headless || true
RUN apt-get -qq -y build-dep hdf5 || true

RUN wget -qO- "https://www.hdfgroup.org/package/source-bzip-2/?wpdmdl=13047&refresh=5c754b75041e11551190901" | tar xjf -


WORKDIR /tmp/hdf5-1.10.4

RUN ./autogen.sh
RUN ./configure --enable-cxx --enable-build-mode=production --prefix=/usr
RUN make
RUN make install

## install H5py separately
RUN apt-get -y -qq install python3-h5py || true
