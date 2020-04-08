FROM ubuntu:18.04

# Set up the basic functionality
ARG DEBIAN_FRONTEND=noninteractive

# Need to install tzdata first for time zones
RUN apt-get update -y \
	&& apt-get install -y tzdata

# Then make sure time zone is set to what it should be
RUN ln -fs /usr/share/zoneinfo/US/Central /etc/localtime \
	&& dpkg-reconfigure --frontend noninteractive tzdata

RUN apt-get update -y \
	&& apt-get install -y gdb language-pack-en-base \
	&& apt-get install -y gcc g++ python2.7 python2.7-dev \
	&& apt-get install -y python2.7-doc python2.7-examples \
	&& apt-get install -y wget vim locate git

RUN apt-get install -y gromacs libgromacs-dev

RUN apt-get install -y libboost1.65-all-dev libboost1.65-dev libboost1.65-doc

RUN apt-get install -y cmake

RUN apt-get install -y cookiecutter
WORKDIR ~/.cookiecutters
RUN git clone https://github.com/mkrapp/cookiecutter-reproducible-science.git

# Set the display to make life easier
ENV DISPLAY=172.50.1.221:0.0
