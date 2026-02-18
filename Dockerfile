#
# A Dockerfile to get RUFUS running 
#

# GCC 4.9 only available up to 16.04
FROM ubuntu:16.04

ARG DEBIAN_FRONTEND=noninteractive 

COPY . /RUFUS

RUN set -ex; \
# Dependencies
	BUILD_DEPS="cmake build-essential g++-4.9 zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev"; \
	apt-get update; \
	apt-get install -y software-properties-common; \
	add-apt-repository ppa:ubuntu-toolchain-r/test; \
	apt-get install -y python wget git bc libgomp1 $BUILD_DEPS; \
# Build RUFUS
	mkdir -p /RUFUS/bin; \
	cd /RUFUS/bin; \
	cmake ../ -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++); \
	make; \
# Build and install htslib 1.11 (provides bgzip and tabix)
	cd /; \
	wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2; \
	tar -xjf htslib-1.11.tar.bz2; \
	cd htslib-1.11; \
	./configure; \
	make; \
	make install; \
	cd /; \
	rm -rf htslib*; \
# Build and install samtools 1.11
	wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2; \
	tar -xjf samtools-1.11.tar.bz2; \
	cd samtools-1.11; \
	./configure; \
	make; \
	make install; \
	cd /; \
	rm -rf samtools*; \
# Verify samtools and htslib installation
	samtools --version; \
	bgzip --version; \
	tabix --version; \
# Install bedtools
	wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary -O /usr/local/bin/bedtools; \
	chmod +x /usr/local/bin/bedtools; \
# Verify bedtools installation
	bedtools --version; \
# Cleanup
	apt-get purge -y --auto-remove $BUILD_DEPS; \
	apt-get clean; \
	rm -rf /var/lib/apt/lists/*; \
	echo done
