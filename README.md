## Multiomics pipeline. PhD project

The proposed pipeline, written mainly in R language is a set of Jupyter templates created to analyse and combine different types of omic data and follows step-by-step protocols designed for each type. The protocols are based on accepted articles and best practices available in the literature.

It is distributed as a docker project that can be locally installed. The provided docker-compose file contains the instructions needed to automatically create a fully working machine with JupyterHub, the pipeline source code and all the necessary libraries and third parties software.

### Instalation 

Only requisite to install the pipeline locally is to have Docker and docker-compose already installed. Follow the instructions in the [project website](https://docs.docker.com/install/). This pipeline has been built using Docker Engine version 19.03.0 and Compose 1.25.5.

Once docker is installed, pull an Ubuntu image from Docker Hub with the command: `sudo docker pull ubuntu:16.04`

Then, clone this proyect to a local directory using `https://github.com/laumadmar/multiomicsPipeline.git` and inspect the content:
- docker-compose.yml: YAML file defining services, networks and volumes
- Dockerfile: text document that contains all the commands a user could call on the command line to assemble an image
- jupyterhub_config.py: configuration file for JupyterHub
- readme.txt: instructions of how to deploy this pipeline
- software: directory containing third parties software

Modify if required the following parameters in docker-compose and/or Jupyter configuration file:
- docker-compose.yml: ports and volumes
- jupyterhub_config.py: c.JupyterHub.base_url = '/jupyter', which contains the base url , and c.JupyterHub.ip=’0.0.0.0’ and c.JupyterHub.port = 8000, the proxy’s main IP address and port. The administrator user is set up in this file using the parameters c.Authenticator.admin_users = {'admin'} and c.JupyterHub.admin_access = True. Where ‘admin’ is the linux user created in the docker-compose.yml file.

### Usage
To build and run locally the pipeline with Docker Compose, follow these steps:
1. From your project directory, start up your application by running `sudo docker-compose up`
Compose builds an image from the instrucctionse specified in the Dockerfile, and starts the services defined. 
2. Enter http://localhost:8000/jupyter in a browser to see the application running. Modify this url conveniently if you have changed the port and the base url in the docker-compose.yml.
3. Press `CTRL+C` to stop the console output. This will also stop the container.
4. From your project directory, keep the app started in the background by running `sudo docker-compose start`

Note you need sudo privileges or create a special group; read more on Docker web site. Few useful docker commands:
- `sudo docker-compose stop` to stop the running container
- `sudo docker-compose ps` to check the status 
- `sudo docker inspect --format='{{.LogPath}}' genomicpipeline_jupyter_1` to get the path to the log file 
- `sudo docker logs genomicpipeline_jupyter_1` to print the log in console
- `sudo docker exec –it genomicpipeline_jupyter_1 bash` to access the running container.

