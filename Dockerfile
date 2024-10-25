FROM rocker/r2u:latest
# 

LABEL author="Angel Angelov <aangeloo@gmail.com>"
LABEL description="Docker image for the NXF-TGS pipeline (merge-rename, faster report, ont wf assemblies)"

RUN apt-get update && apt-get install -y \
    pandoc nano ksh procps libxt-dev libssl-dev libxml2-dev libfontconfig1-dev parallel python3-pip pigz ed
# libxt-dev is required to solve the segfault error caused by cairoVersion() in R

# RUN curl -s https://get.nextflow.io | bash && chmod +x nextflow && mv nextflow /usr/local/bin/

# setup faster and fastkmers for linux
RUN wget -P bin https://github.com/angelovangel/faster/releases/download/v0.2.1/x86_64-linux-faster && \
    mv bin/x86_64-linux-faster bin/faster && \
    chmod 755 bin/faster

RUN wget -P bin https://github.com/angelovangel/faster2/releases/download/v0.3.0/faster2 && \
    chmod 755 bin/faster2

RUN wget -P bin https://github.com/angelovangel/fasterplot/releases/download/v0.1.0/fasterplot && \
    chmod 755 bin/fasterplot

RUN wget -P bin https://github.com/angelovangel/fastkmers/releases/download/v0.1.3/fastkmers && \
    chmod 755 bin/fastkmers

# seqkit
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz && \
    tar -zxvf seqkit_linux_amd64.tar.gz && \
    cp seqkit bin/

# samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
    tar -jxf samtools-1.20.tar.bz2 && \
    rm samtools-1.20.tar.bz2 && \
    cd samtools-1.20 && \
    ./configure --prefix $(pwd) && \
    make

ENV PATH=${PATH}:/samtools-1.20

# minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 && \
    tar -jxf minimap2-2.28_x64-linux.tar.bz2 && \
    rm minimap2-2.28_x64-linux.tar.bz2 
    
ENV PATH=${PATH}:/minimap2-2.28_x64-linux

# perbase
RUN wget -P bin https://github.com/sstadick/perbase/releases/download/v0.9.0/perbase-linux-amd64 && \
    mv bin/perbase-linux-amd64 bin/perbase && \
    chmod 755 bin/perbase

# igv-reports
RUN pip install igv-reports

RUN install2.r \
    'R.utils' \
    stringr \
    readr \
    readxl \
    writexl \
    knitr \
    DT \
    kableExtra \
    dplyr \
    sparkline \
    htmlwidgets \
    jsonlite \
    #parallel \
    parallelMap \
    optparse \
    rmarkdown \
    funr \
    && rm -rf /tmp/downloaded_packages
