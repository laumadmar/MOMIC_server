FROM ubuntu:18.04

ENV HTTP_PROXY "http://192.168.20.101:3128"
ENV HTTPS_PROXY "http://192.168.20.101:3128"
ENV DEBIAN_FRONTEND=noninteractive

# Updating Ubuntu packages
RUN apt-get update && yes|apt-get upgrade

# Adding vim and other necessary libraries (texlive-xxx for dowloading notebooks as pdf)
RUN apt-get update --fix-missing && apt-get install -y --no-install-recommends curl vim wget git tabix software-properties-common texlive-latex-extra lmodern vcftools zlib1g-dev liblzma-dev libbz2-dev gfortran libxml2-dev libncurses5-dev libfontconfig1-dev libcairo2-dev dirmngr gpg-agent

# Enable the CRAN repository and add the CRAN GPG key for installing R 4.0. For older versions this is not needed.
#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN echo "deb http://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/" > /etc/apt/sources.list.d/cran.list

# Install R 4.0
RUN apt-get update && apt-get install -y r-base r-base-dev r-recommended

# Installing python, pip and node.js (use with pip3 or python3 -m pip)
#RUN apt-get install -y python3 python3-pip 
RUN apt-get install -y python3-pip
RUN curl -sL https://deb.nodesource.com/setup_12.x |  bash -
RUN apt-get install -y nodejs

# Installing jupyterhub
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip --version
RUN python3 -m pip install jupyterhub
RUN npm install -g configurable-http-proxy
RUN python3 -m pip install notebook  

# Install jupyterLab, git integration, jupyter extensions and extensions_configurator to be able to see the extensions in the notebook
RUN python3 -m pip install --upgrade jupyterlab jupyterlab-git==0.30.0b2
#RUN jupyter lab build
RUN jupyter labextension install @jupyterlab/git
RUN jupyter labextension install @jupyterlab/toc

RUN python3 -m pip install jupyter_contrib_nbextensions 
RUN jupyter contrib nbextension install --system

RUN jupyter nbextension enable execute_time/ExecuteTime 
RUN jupyter nbextension enable toc2/main
RUN jupyter nbextension enable varInspector/main

# R kernel
RUN apt-get install -y libzmq3-dev libcurl4-openssl-dev libssl-dev
RUN R -e "install.packages(c('repr', 'IRdisplay', 'IRkernel'), type = 'source',repos='http://cran.rstudio.com/')"
RUN R -e "IRkernel::installspec(user = FALSE)"

# Install rpy2
RUN apt-get install libffi-dev
RUN python3 -m pip install rpy2
RUN python3 -m pip install numpy
RUN python3 -m pip install pandas

# copying software folder to install third-party software
COPY software /tmp

# Installing R packages
#RUN Rscript /tmp/restore_R_packages.R
#RUN R -e "install.packages('BiocManager')"
#RUN R -e "BiocManager::install(c('BiocGenerics','Biobase','affy','gcrma','limma','hgu133plus2cdf','hgu133plus2.db','hgu133a.db','AnnotationDbi','DESeq2','biomaRt','lumi','GEOquery','preprocessCore','sva','edgeR','impute','DEqMS'), update = TRUE, ask = FALSE)"

#RUN R -e "install.packages(c('RSQLite','ploty','qqman','gridExtra','grid','gplots','NMF','ggplot2','knitr','WebGestaltR','pheatmap','calibrate','GOplot','RobustRankAggreg'))"
#RUN R -e "library(remotes); remotes::install_github('metaOmics/MetaDE')"

# Installing Java
RUN mkdir /usr/java && \
tar -xzvf /tmp/jre-8u121-linux-x64.tar.gz -C /usr/java

# Installing EIG software
RUN tar -xvzf /tmp/EIG-6.1.4.tar.gz -C /usr/lib

#Installing PLINK
RUN unzip /tmp/plink1.9.zip -d /usr/lib/plink1.9 && \
dpkg -i /tmp/pandoc-1.19.2-1-amd64.deb

# Installing magma
RUN unzip /tmp/magma_v1.06b.zip -d /usr/lib/magma && \
chmod +x /usr/lib/magma/magma

# Installing METAL
RUN tar -xvzf /tmp/Linux-metal.tar.gz -C /usr/lib

# Install FASTQC and multiqc
RUN unzip -q fastqc_v0.11.9.zip -d /usr/lib/ && \
chmod +x /usr/lib/FastQC/fastqc
RUN python3.6 -m pip install multiqc 

# Install STAR aligner
RUN tar -xzf /tmp/2.7.7a.tar.gz -C /usr/lib
# RUN cd /usr/lib/STAR-2.7.7a/source && \
# RUN make STAR && cd

# Install libpng12-0 needed for liftOver
RUN dpkg -i /tmp/libpng12-0_1.2.54-1ubuntu1.1_amd64.deb

# setting up path variable
ENV PATH="/usr/lib/EIG-6.1.4/bin:/usr/java/jre1.8.0_121/bin:/usr/lib/plink1.9:/usr/lib/magma:/usr/lib/generic-metal:/usr/lib/FastQC/fastqc:/usr/lib/STAR-2.7.7a/bin/Linux_x86_64:${PATH}"

# create directory where to keep jupyter config files
RUN mkdir /opt/jupyterhub

# add jupyterhub_config.py  
ADD jupyterhub_config.py /opt/jupyterhub/jupyterhub_config.py

# create directory where to keep pipeline files
RUN mkdir -p /mnt/data

# Add momic user
RUN useradd -u 1005 --create-home --shell /bin/bash momic
RUN echo 'momic:m0m1c' | chpasswd
RUN chown momic:momic /home/momic

# start jupyterhub as root
USER root
WORKDIR /opt/jupyterhub
ENTRYPOINT ["jupyterhub"]
