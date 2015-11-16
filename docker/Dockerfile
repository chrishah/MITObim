FROM ubuntu:14.04

MAINTAINER Christoph Hahn <chrisi.hahni@gmail.com>

RUN apt-get update
RUN apt-get upgrade -y

# Not essential, but wise to set the lang (taken from szitenberg/reprophylo)
RUN apt-get install -y language-pack-en
ENV LANGUAGE en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8

RUN locale-gen en_US.UTF-8
RUN dpkg-reconfigure locales

# some necessary dependencies, developer tools
RUN apt-get install -y -q build-essential make gcc zlib1g-dev git python python-dev python-pip
RUN apt-get install git

# add scripts, executables, testdata
RUN mkdir /home/src
ADD external_software /home/src/
RUN mkdir /home/src/scripts
ADD scripts/MITObim_1.8.pl /home/src/scripts/
ADD scripts/MITObim_1.6.pl /home/src/scripts/
ADD scripts/interleave-fastqgz-MITOBIM.py /home/src/scripts/
RUN chmod a+x /home/src/scripts/MITObim_1.8.pl
RUN chmod a+x /home/src/scripts/MITObim_1.6.pl
RUN chmod a+x /home/src/scripts/interleave-fastqgz-MITOBIM.py
ADD testdata/testdata1.tgz /home/testdata/
RUN mkdir /home/temp

RUN ln -s /home/src/mira_4.0.2/mira* /usr/bin
ENV PATH /home/src/scripts:$PATH
# RUN ln -s /home/src/scripts/MITObim_1.8.pl /usr/bin

VOLUME /home/data
WORKDIR /home
RUN chmod -R a+rw /home

