FROM rocker/r2u:latest

LABEL author="Angel Angelov <aangeloo@gmail.com>"
LABEL description="Docker image for the NXF-TGS pipeline (merge-rename, faster report, ont wf assemblies)"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        pandoc nano ksh procps libxt-dev libssl-dev libxml2-dev libfontconfig1-dev \
        parallel pipx git curl cmake pigz ed build-essential wget ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Install Rust and build binaries, then clean up
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y && \
    export PATH="/root/.cargo/bin:${PATH}" && \
    mkdir -p /bin && \
    for repo in fastkmers faster faster2 fasterplot; do \
      git clone --depth 1 https://github.com/angelovangel/$repo.git && \
      cd $repo && \
      cargo build --release && \
      cp target/release/$repo /bin/ || true && \
      cd .. && rm -rf $repo; \
    done && \
    rm -rf /root/.cargo/registry /root/.cargo/git

# seqkit
RUN wget -q https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz && \
    tar -zxf seqkit_linux_amd64.tar.gz && \
    mv seqkit /bin/ && \
    rm seqkit_linux_amd64.tar.gz

# samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
    tar -jxf samtools-1.20.tar.bz2 && \
    cd samtools-1.20 && \
    ./configure --prefix=$(pwd) && make && \
    cp samtools ../bin/ && \
    cd .. && rm -rf samtools-1.20 samtools-1.20.tar.bz2

# minimap2 (build from source)
RUN git clone --depth 1 https://github.com/lh3/minimap2.git && \
    cd minimap2 && make && \
    cp minimap2 ../bin/ && \
    cd .. && rm -rf minimap2

# perbase
RUN wget -q -P /bin https://github.com/sstadick/perbase/releases/download/v0.9.0/perbase-linux-amd64 && \
    mv /bin/perbase-linux-amd64 /bin/perbase && \
    chmod 755 /bin/perbase

# igv-reports
ENV PATH="/root/.local/bin:/bin:${PATH}"
RUN pipx ensurepath && pipx install igv-reports

RUN install2.r \
    'R.utils' \
    stringr \
    readr \
    readxl \
    knitr \
    bslib \
    bsicons \
    shiny \
    scales \
    reactable \
    DT \
    kableExtra \
    dplyr \
    sparkline \
    htmlwidgets \
    jsonlite \
    parallelMap \
    optparse \
    rmarkdown \
    funr \
    vroom \
    && rm -rf /tmp/downloaded_packages
