# Installation

Before installation, look into the [[index#PC requirements|PC requirements]] if the system is suitable for running the dashboard tool.
## Set up the server
MetaPepView runs on top of the dash framework in python. To run the dashboard, the server needs to be started locally on the PC. The simples approach is to start a [Docker](https://www.docker.com/) container from the DockerFile provided in the Github repository. For this, docker (or another OCI compliant manager) has to be installed on the system. Alternatively, the server can be started by setting up the python environment and run the entrypoint file.
### Docker setup.

First, a [Docker image](https://docs.docker.com/get-started/docker-concepts/building-images/build-tag-and-publish-an-image/) needs to be built:

!!! note
    Make sure the docker daemon is running during the command execution (start docker desktop if installed).

```Bash
$ docker build -t Metapepview . # directory of Dockerfile
```

Once the image is successfully built, start a container by running the following command:

```Bash
$ docker run -p 8050:8050 Metapepview
```

`-p` publishes port 8050 from the container to port 8050 in the host pc. This allows us to reach the dashboard in the container. We can access the dashboard by typing in `http://localhost:8050` in the web browser. Note that the dashboard is only accessible from the host PC as long as the port is closed in the firewall.

### Python environment setup

The repository provides a requirements file to install all required python packages with [pip](https://packaging.python.org/en/latest/tutorials/installing-packages/#requirements-files). To execute the tool directly on the host pc, a recent version of python needs to be installed (at least 3.11). Then, The dashboard can be run directly on the host pc by executing:

```Bash
$ pip install -r ./requirements.txt
$ python ./index.py
```