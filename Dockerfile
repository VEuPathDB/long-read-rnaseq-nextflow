FROM biocontainers/samtools:v1.9-4-deb_cv1
FROM python:3.7

# Set environment variables
ENV talon_version=5.0
ENV bedtools_version=2.30.0
ENV TranscriptClean_version=2.0.3
ENV DEBIAN_FRONTEND=noninteractive


USER root
# Download dependencies
RUN apt-get update && apt-get install -y build-essential wget unzip default-jre tabix procps && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*
RUN pip install --upgrade pip
RUN apt-get update && apt-get install --no-install-recommends -y libncurses5-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev zlib1g-dev libssl-dev gcc make perl bzip2 gnuplot ca-certificates gawk && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*
WORKDIR /usr/local/

# download software
RUN wget  https://github.com/mortazavilab/TALON/archive/refs/tags/v${talon_version}.zip 
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${bedtools_version}/bedtools-${bedtools_version}.tar.gz  
RUN wget https://github.com/mortazavilab/TranscriptClean/archive/refs/tags/v${TranscriptClean_version}.zip 



# unpack
RUN unzip v${talon_version}.zip 
RUN tar -zxvf bedtools-${bedtools_version}.tar.gz 
RUN unzip v${TranscriptClean_version}.zip 


WORKDIR /usr/local/TALON-${talon_version}
RUN pip install cython 
RUN pip install "setuptools<58.0.0" 
RUN pip install . 

WORKDIR /usr/local/bedtools2
RUN make 
RUN ln -s /usr/local/bedtools2/bin/bamToBed /usr/local/bin/bamToBed 
RUN ln -s /usr/local/bedtools2/bin/bedtools /usr/local/bin/bedtools


RUN pip install pyfaidx 
RUN pip install pyranges
RUN pip install pyfasta
RUN pip install pandas




WORKDIR /usr/local/TranscriptClean-${TranscriptClean_version}
RUN cp -r ./* /usr/local/bin/
COPY ./bin/* /usr/local/bin/

WORKDIR /work

