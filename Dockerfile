FROM rocker/r2u:latest
# 

LABEL author="Angel Angelov <aangeloo@gmail.com>"
LABEL description="Docker image for the NXF-TGS pipeline (merge-rename, faster report, ont wf assemblies)"

RUN apt-get update && apt-get install -y \
    pandoc nano ksh procps libxt-dev libssl-dev libxml2-dev libfontconfig1-dev parallel
# libxt-dev is required to solve the segfault error caused by cairoVersion() in R

# RUN curl -s https://get.nextflow.io | bash && chmod +x nextflow && mv nextflow /usr/local/bin/

# setup faster and fastkmers for linux
RUN wget -P bin https://github.com/angelovangel/faster/releases/download/v0.2.1/x86_64-linux-faster && \
    mv bin/x86_64-linux-faster bin/faster && \
    chmod 755 bin/faster

RUN wget -P bin https://github.com/angelovangel/faster2/releases/download/v0.3.0/faster2 && \
    chmod 755 bin/faster2

RUN wget -P bin https://github.com/angelovangel/fastkmers/releases/download/v0.1.3/fastkmers && \
    chmod 755 bin/fastkmers

RUN install2.r \
    'R.utils' \
    stringr \
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
