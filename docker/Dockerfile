FROM ubuntu:16.04

MAINTAINER Christoph Hahn <christoph.hahn@uni-graz.at>

RUN apt-get update && apt-get upgrade -y

# Not essential, but wise to set the lang (taken from szitenberg/reprophylo)
RUN apt-get install -y language-pack-en
ENV LANGUAGE en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8

RUN locale-gen en_US.UTF-8 && dpkg-reconfigure locales

# some necessary dependencies, developer tools
# set up the latest GNU parallel
RUN apt-get install -y -q build-essential git wget curl python python-dev python-pip && \
wget https://ftpmirror.gnu.org/parallel/parallel-20180822.tar.bz2 && \
bzip2 -dc parallel-20180822.tar.bz2 | tar xvf - && \
cd parallel-20180822 && \
./configure && make && make install && cd .. && rm -rf parallel-20180822
 
#&& \
#(wget -O - pi.dk/3 || curl pi.dk/3/ || fetch -o - http://pi.dk/3) | bash

# add scripts, executables, testdata (the target directories will be made automatically if they don't exist)
ADD external_software /home/src/
ADD scripts /home/src/scripts/
ADD testdata/testdata1.tgz /home/testdata/

RUN chmod -R a+x /home/src/scripts && ln -s /home/src/mira_4.0.2/mira* /usr/bin
ENV PATH /home/src/scripts:$PATH

VOLUME /home/data
WORKDIR /home
RUN chmod -R a+rw /home

