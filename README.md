## MOMIC: A Multi-Omics Pipeline for data analysis, integration and interpretation

MOMIC offers a complete analysis environment for analysing and integrating multi-omics data in a single, easy-to-use platform.

MOMIC currently compiles protocols for whole genome SNP data (GWAS), mRNA expression (both from arrays and from RNAseq experiments) and protein data. Along with enrichment analysis and methods for combining heterogeneous data at different molecular levels. The proposed protocols are developed as Jupyter notebooks that guide the user through the tasks of pre-processing and transforming the data and performing the actual analysis, allowing the user to modify any piece of code needed along the process to adequate it to each project.

It is distributed as a docker project and a collection of Jupyter notebooks. The docker-compose file provided in this repository contains the instructions needed to automatically create a fully working machine with JupyterHub, the pipeline source code and all the necessary libraries and third parties software.

### Instalation 

Only requisite to install the pipeline locally is to have Docker and docker-compose already installed. Follow the instructions in the [project website](https://docs.docker.com/install/). This pipeline has been built using Docker Engine version 19.03.0 and Compose 1.25.5. Notice you need sudo privileges or a special group for running docker commands; read more on the Docker web site. 

Once docker is installed, pull MOMIC image from Docker Hub: `docker pull laumadmarq/momic:latest`. You can ensure that the image is installed by using `docker images`.

Then, clone this project to a local directory using `https://github.com/laumadmar/MOMIC_server.git` and inspect the content:
- docker-compose.yml: YAML file defining services, networks and volumes
- Dockerfile: text file that contains the instructions to build the service
- jupyterhub_config.py: configuration file for JupyterHub
- README.md: instructions for deploying this pipeline
- software: directory containing third-party software

Modify if required the following parameters in docker-compose and/or Jupyter configuration file:
- docker-compose.yml: ports and volumes
- jupyterhub_config.py: c.JupyterHub.base_url = '/jupyter', which contains the base url , and c.JupyterHub.ip=0.0.0.0â and c.JupyterHub.port = 8000, the proxy's main IP address and port. The administrator user is set up in this file using the parameters c.Authenticator.admin_users = {'admin'} and c.JupyterHub.admin_access = True. Where admin is the linux user created in the docker-compose.yml file.

It is advisable to check out the [user manual](https://laumadmar.github.io/MOMIC_server) to get a better understanding of the instalation process and customisation options.

#### Docker instructions
Once you have MOMIC docker image and server, build and run the service locally following these steps from the command line:
1. From your project directory, start up your application by running `sudo docker-compose up`
Compose builds an image from the instrucctionse specified in the Dockerfile, and starts the services defined. 
2. Once the server is up, press `CTRL+c` to stop the console output. This will also stop the container.
3. From your project directory, keep the app started in the background by running `sudo docker-compose start`
5. Enter http://localhost:8000/jupyter in a browser to see the application running. Modify this url conveniently if you have changed the port and the base url in the docker-compose.yml. Log in with user momic, pass m0m1c.

An alternative is to create your container from the original instructions, instead of using momic docker image, which can be fully customised. After clonning MOMIC_server from the repository, rename the file `Dockerfile.steps` to `Dockerfile`.
1. run `docker pull ubuntu:18.04` to download the ubuntu docker image
2. run `docker-compose up` and once the server is up, press `CTRL+c` to stop the console output
3. run `docker-compose start` to keep the service running in the background
4. ssh into the container executing the access script (`./access`) and type `nohup Rscript /tmp/install_specific_libraries.R &` to install the required R packages
5. access the tool at http://localhost:8000/jupyter and log in
6. clone MOMIC Notebooks repo from jupyter or from the terminal after ssh into the container. For the former, go to the git tab located in the left menu, click on the button `Clone a Repository` and provide the repo url. For the latter, ssh into the container, cd into momic home directory and type  `git clone https://github.com/laumadmar/MOMIC_notebooks.git`

#### Access and log scripts
Create a bash script to access the machine and another one to get the logs:

`docker exec -it momic_server_web_1 bash`

`docker logs momic_server_web_1`

Change momic_server_web_1 by the name of the service if it differs. You can get if from `docker-compose ps`.

Note you need sudo privileges or create a special group; read more on Docker web site. Few useful docker commands:
- `docker-compose stop` to stop the running container
- `docker-compose ps` to check the status
- `docker inspect --format='{{.LogPath}}' momic_server_1` to get the path to the log file
- `docker logs momic_server_web_1` to print the log in console
- `docker exec -it momic_server_web_1 bash` to access the running container.

