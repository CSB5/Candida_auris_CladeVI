FROM ubuntu:22.04
RUN apt-get update && apt-get upgrade -y
RUN apt-get install libssl-dev
RUN apt-get install -y -q build-essential curl wget unzip

# install chopper
# Chopper require GLIBC_2.32 or GLIBC_2.33 or GLIBC_2.34 which are in a newer kernel
RUN wget https://github.com/wdecoster/chopper/releases/download/v0.5.0/chopper-linux.zip \
    && unzip chopper-linux.zip \
	&& chmod a+x chopper \
	&& mv chopper /usr/local/bin
