FROM bioconductor/bioconductor_docker:RELEASE_3_19
WORKDIR /app #create directory
RUN apt-get update && apt-get install -y \ 
      libhdf5-dev \          
      libcurl4-openssl-dev \ 
      libssl-dev \           
      libxml2-dev \         
      libgit2-dev \         
      && rm -rf /var/lib/apt/lists/*  
COPY install.packages.R .
RUN Rscript install.packages.R 
COPY scripts/ scripts/ 
USER rstudio

#The bioconductor_docker image is built on top of Ubuntu
# libhdf5-dev  to build hdf5r — lets Seurat read .h5 files
# libcurl4-openssl-dev \ # to make internet connections — used when downloading packages
# libssl-dev \ # for secure connections (https) — same reason
# libxml2-dev \  # to parse XML — needed by some Bioconductor packages
# libgit2-dev \ # to interact with git — needed by remotes to install from GitHub
# && rm -rf /var/lib/apt/lists/*  # clean up download cache to keep image size small
#COPY <from /scripts on my machine> to </scripts in container>

