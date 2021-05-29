## MOMIC: A Multi-Omics Pipeline for data analysis, integration and interpretation

MOMIC offers a complete analysis environment for analysing and integrating multi-omics data in a single, easy-to-use platform.

MOMIC currently compiles protocols for whole genome SNP data (GWAS), mRNA expression (both from arrays and from RNAseq experiments) and protein data. The proposed protocols are developed as Jupyter notebooks that guide the user through the tasks of pre-processing and transforming the data and performing the actual analysis, allowing the user to modify any piece of code needed along the process to adequate it to each project.

It is distributed as a docker project and a collection of Jupyter notebooks. The docker-compose file provided in this repository contains the instructions needed to automatically create a fully working machine with JupyterHub, the pipeline source code and all the necessary libraries and third parties software. Once you have your local MOMIC server up and running, get the notebooks following the instrucctions on [this repo](https://github.com/laumadmar/MOMIC_notebooks.git) or check out the [user manual](https://laumadmar.github.io/MOMIC_server).

### Instalation 

Only requisite to install the pipeline locally is to have Docker and docker-compose already installed. Follow the instructions in the [project website](https://docs.docker.com/install/). This pipeline has been built using Docker Engine version 19.03.0 and Compose 1.25.5.

Once docker is installed, pull an Ubuntu image from Docker Hub with the command: `sudo docker pull ubuntu:18.04`

Then, clone this proyect to a local directory using `https://github.com/laumadmar/MOMIC_server.git` and inspect the content:
- docker-compose.yml: YAML file defining services, networks and volumes
- Dockerfile: text document that contains all the commands a user could call on the command line to assemble an image
- jupyterhub_config.py: configuration file for JupyterHub
- README.md: instructionsfor deploying this pipeline
- software: directory containing third-party software

Modify if required the following parameters in docker-compose and/or Jupyter configuration file:
- docker-compose.yml: ports and volumes
- jupyterhub_config.py: c.JupyterHub.base_url = '/jupyter', which contains the base url , and c.JupyterHub.ip=0.0.0.0â and c.JupyterHub.port = 8000, the proxy's main IP address and port. The administrator user is set up in this file using the parameters c.Authenticator.admin_users = {'admin'} and c.JupyterHub.admin_access = True. Where admin is the linux user created in the docker-compose.yml file.

#### Docker instructions
To build and run locally the pipeline with Docker Compose, follow these steps:
1. From your project directory, start up your application by running `sudo docker-compose up`
Compose builds an image from the instrucctionse specified in the Dockerfile, and starts the services defined. 
2. Enter http://localhost:8000/jupyter in a browser to see the application running. Modify this url conveniently if you have changed the port and the base url in the docker-compose.yml.
3. Press `CTRL+C` to stop the console output. This will also stop the container.
4. From your project directory, keep the app started in the background by running `sudo docker-compose start`

**Important Note:** The instalation of the necesary R libraries is commented out in the Dockerfile as it takes a lot to build. You can either uncomment these lines or access the container once started and run this on the terminal `nohup Rscript /tmp/install_specific_libraries.R &`. 

#### Access and log scripts
Create a bash script to access the machine and another one to get the logs:

sudo docker exec -it jupyter_config_web_1 bash

sudo docker logs jupyter_config_web_1

Change jupyter_config_web_1 by the name of the service. You can get if from sudo docker-compose ps.

Note you need sudo privileges or create a special group; read more on Docker web site. Few useful docker commands:
- `sudo docker-compose stop` to stop the running container
- `sudo docker-compose ps` to check the status 
- `sudo docker inspect --format='{{.LogPath}}' genomicpipeline_jupyter_1` to get the path to the log file 
- `sudo docker logs genomicpipeline_jupyter_1` to print the log in console
- `sudo docker exec -it genomicpipeline_jupyter_1 bash` to access the running container.

