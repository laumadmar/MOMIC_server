FROM laumadmarq/momic:latest
  
# Updating Ubuntu packages
RUN apt-get update && yes|apt-get upgrade

# start jupyterhub as root
USER root
WORKDIR /opt/jupyterhub
ENTRYPOINT ["jupyterhub"]
