FROM ubuntu:latest
LABEL dockerfile.version="1"
LABEL software="bcftools, samtools, PLINK, R"
#LABEL license="https://github.com/.../LICENSE"
LABEL creator="Dafni Michalettou"
LABEL creator.email="dafni.michalettou@newcastle.ac.uk"

ARG samtoolsVer="1.21"
ARG bcftoolsVer="1.21"
#ARG minicondaVer="311_25.1.1-2"
#If using plink 1.9 ARGs must change a bit
ARG PLINK_v="2"
ARG PLINK_date="20250420"
ARG R_VERSION="4.4.3"

#   Install basics [tabix includes bgzip]
#ARG DEBIAN_FRONTEND=noninteractive
#ENV TZ=Europe/London

RUN apt-get update && apt-get install --no-install-recommends -y \
    wget build-essential zlib1g-dev liblzma-dev libbz2-dev libxau-dev libgsl-dev libpcre2-dev && \
    apt-get install -y curl make automake autoconf gcc gfortran \
    bzip2 unzip tabix \
    perl ca-certificates \
    libncurses5-dev libncursesw5-dev libcurl4-gnutls-dev libssl-dev libperl-dev libgsl0-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get -y install tzdata &&\
    apt-get -y clean && \
    rm -rf /var/lib/apt/lists/* && apt-get autoclean
#   python3-pip
#   apt install build-essential software-properties-common libssl-dev libffi-dev python3-dev libgdbm-dev libc6-dev libbz2-dev libsqlite3-dev tk-dev libffi-dev zlib1g-dev -y

##   Install HTSLIB
##    cd /usr/bin
#RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
#    tar -vxjf htslib-1.9.tar.bz2 && \
#    cd htslib-1.9 && \
#    make

#   Install SAMTOOLS
RUN wget https://github.com/samtools/samtools/releases/download/${samtoolsVer}/samtools-${samtoolsVer}.tar.bz2 && \
    tar -vxjf samtools-${samtoolsVer}.tar.bz2 && \
    rm samtools-${samtoolsVer}.tar.bz2 &&\
    cd samtools-${samtoolsVer} &&\
    make &&\
    make install &&\
    cd ..

#   Install BCFTools
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftoolsVer}/bcftools-${bcftoolsVer}.tar.bz2 && \
    tar -vxjf bcftools-${bcftoolsVer}.tar.bz2 && \
    rm bcftools-${bcftoolsVer}.tar.bz2 && \
    cd bcftools-${bcftoolsVer} &&\
    make &&\
    make install &&\
    cd ..

#   Export To Path And Refresh
#   export PATH="$PATH:/usr/bin/bcftools-1.9"
#   export PATH="$PATH:/usr/bin/samtools-1.9"
#   export PATH="$PATH:/usr/bin/htslib-1.9"
#   source ~/.profile

##  Install miniconda [to install PLINK]
#RUN mkdir -p /opt/conda 
#RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py${minicondaVer}-Linux-x86_64.sh -O /opt/conda/miniconda.sh && \
#    bash /opt/conda/miniconda.sh -b -p /opt/conda

# Install your environment.yaml deps into base env 
#COPY environment.yaml /tmp 
#RUN . /opt/miniconda/bin/activate && conda env update --name base --file /tmp/environment.yaml 

#   Install PLINK
RUN wget https://s3.amazonaws.com/plink${PLINK_v}-assets/alpha5/plink${PLINK_v}_linux_x86_64_${PLINK_date}.zip && \
    unzip plink${PLINK_v}_linux_x86_64_${PLINK_date}.zip -d /usr/bin/ && \
    rm -v plink${PLINK_v}_linux_x86_64_${PLINK_date}.zip

#     wget https://www.cog-genomics.org/static/bin/plink$PLINK_VERSION/plink_linux_x86_64.zip && \
#     unzip plink_linux_x86_64.zip && \
#     rm plink_linux_x86_64.zip && \

#   Install R
RUN wget -c https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz && \
    tar -xf R-${R_VERSION}.tar.gz && \
    cd R-${R_VERSION} && \
    ./configure --with-readline=no --with-x=no && \
#    make -j$(nproc) && \
    make &&\
    make install && \
    cd .. && \
    rm -rf R-${R_VERSION} R-${R_VERSION}.tar.gz

#RUN apt-get update && \
#    apt-get install -y --no-install-recommends software-properties-common dirmngr
#RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
#RUN sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
#RUN apt-get update && \
#    apt-get install -y r-base r-base-core r-recommended r-base-dev
