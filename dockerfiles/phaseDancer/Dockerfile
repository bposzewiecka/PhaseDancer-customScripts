FROM ubuntu:focal

RUN apt-get update -y
RUN apt-get upgrade -y 

RUN apt-get install git -y
RUN apt-get install make -y
RUN apt-get install build-essential -y
RUN apt-get install libz-dev -y
RUN apt-get install autoconf autoconf-archive -y

RUN apt-get install python3-pip -y
RUN apt-get install zip -y
RUN apt-get install automake gcc libbz2-dev liblzma-dev -y

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata

RUN apt-get install software-properties-common -y
RUN add-apt-repository ppa:deadsnakes/ppa -y
RUN apt-get update -y
RUN apt-get install python3.8 -y
RUN apt-get install python-is-python3 -y
RUN apt-get install gawk -y

RUN apt-get install samtools -y

RUN apt-get install expect -y

RUN mkdir phaseDancer_bin

###################################### MINIMAP ######################################

WORKDIR /phaseDancer_bin
RUN git clone https://github.com/lh3/minimap2
WORKDIR /phaseDancer_bin/minimap2 
RUN make
ENV PATH="/phaseDancer_bin/minimap2:${PATH}"

###################################### WTDBG2 #######################################

WORKDIR /phaseDancer_bin
RUN git clone https://github.com/ruanjue/wtdbg2
WORKDIR /phaseDancer_bin/wtdbg2
RUN make
ENV PATH="/phaseDancer_bin/wtdbg2:${PATH}"

#################################### SNAKEMAKE ######################################

RUN pip install snakemake

#################################### PHASEDANCER ####################################

RUN pip install numpy
RUN pip install pysam
RUN pip install pyfaidx
RUN pip install biopython
RUN pip install edlib
RUN pip install sklearn
RUN pip install matplotlib
RUN pip install scipy

WORKDIR /phaseDancer_bin
RUN git clone https://github.com/bposzewiecka/phaseDancer.git
ENV PHASEDANCER_DIR=/phaseDancer_bin/phaseDancer

RUN mkdir -p /phaseDancerData/data
ENV PHASEDANCER_DATA_DIR=/phaseDancerData

WORKDIR /phaseDancer_bin/phaseDancer

ENTRYPOINT ["bash"]
