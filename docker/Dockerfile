FROM ubuntu:14.04

MAINTAINER Christoph Hahn <chrisi.hahni@gmail.com>

RUN apt-get update
RUN apt-get upgrade -y

# Not essential, but wise to set the lang
RUN apt-get install -y language-pack-en
ENV LANGUAGE en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8

RUN locale-gen en_US.UTF-8
RUN dpkg-reconfigure locales

# Python binary dependencies, developer tools
RUN apt-get install -y -q build-essential make gcc zlib1g-dev git python python-dev python-pip

RUN apt-get install git

RUN mkdir /home/src
ADD external_software /home/src/
RUN mkdir /home/src/scripts
ADD scripts/MITObim_1.8.pl /home/src/scripts/
RUN chmod a+x /home/src/scripts/MITObim_1.8.pl
ADD scripts/MITObim_1.6.pl /home/src/scripts/
RUN chmod a+x /home/src/scripts/MITObim_1.6.pl
ADD testdata/testdata1.tgz /home/testdata/

RUN ln -s /home/src/mira_4.0.2/mira* /usr/bin
RUN ln -s /home/src/scripts/MITObim_1.8.pl /usr/bin

VOLUME /home/data
WORKDIR /home
RUN chmod -R a+rw /home

