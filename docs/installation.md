# Installation

Before installation, look into the [[index#PC requirements|PC requirements]] if the system is suitable for running the dashboard tool.
## Set up the server
MetaPepView runs on top of the dash framework in python. To run the dashboard, the server needs to be started locally on the PC. The simplest approach is to start a [Docker](https://www.docker.com/) container from the DockerFile provided in the Github repository. For this, Docker (or another OCI compliant manager) has to be installed on the system. Alternatively, MetaPepView can be installed with `pipx`.

### Run in Docker.

First, a [Docker image](https://docs.docker.com/get-started/docker-concepts/building-images/build-tag-and-publish-an-image/) needs to be built:

!!! note
    Make sure the Docker daemon is running during the command execution (start Docker desktop if installed).

```Bash
$ docker build -t Metapepview . # directory of Dockerfile
```

Once the image is successfully built, start a container by running the following command:

```Bash
$ docker run -p 8050:8050 Metapepview
```

the options`-p` publishes port 8050 from the container to port 8050 in the host pc. This allows us to reach the dashboard in the container. We can access the dashboard in the web browser at  `http://localhost:8050`. Note that the dashboard is only accessible from the host PC as long as the port is closed in the firewall.

### Install with `pipx` (or `pip`)

The full MetaPepView application can be installed on any PC system with `pip`. This installation provides the MetaPepView dashboard, as well as a command line utility for construction of [[experiment-evaluation|experimental reference datasets]]. It is highly recommended to install the application with [pipx](https://pipx.pypa.io/stable/). This provides a separate environment for the tool to run without cluttering existing python environments and ensures the executables are exposed to the command line.

To install MetaPepView, first, install [python](https://www.python.org/).
Once installed, [install pipx](https://pipx.pypa.io/stable/installation/).

#### Installation with `pipx`

In Ubuntu:
```Bash
sudo apt install pipx
pipx ensurepath
```

In Fedora:
```Bash
sudo dnf install pipx
pipx ensurepath
```

In windows, `pipx` can be installed via `pip`:
```Powershell
# If you installed python using Microsoft Store, replace `py` with `python3` in the next line.
py -m pip install --user pipx
python -m pipx ensurepath
```

Once pipx is installed, MetaPepView can be installed:
```Bash
pipx install metapepview
```

Once installed, the MetaPepView server can be started by running the application in the command line:
```Bash
metapepview
```

The dashboard can be accessed in the web browser.

To uninstall MetaPepView:
```Bash
pipx uninstall metapepview
```

#### Installation with `pip` (Not recommended)

When python is installed, MetaPepView may be installed directly with `pip`. Installing this way
will add MetaPepView as package to the global python environment. In addition, it is not certain
that the application executables are exposed on the command line. To provide some isolation, it 
is recommended to install MetaPepView in a [virtual environment](https://docs.python.org/3/library/venv.html)
or in a [conda environment](https://docs.conda.io/en/latest/)

To install MetaPepView with pip:
```Bash
pip install metapepview
```

If the python environment is exposed to the operating system PATH, the MetaPepView server can be started by running the application in the command line:
```Bash
metapepview
```

To uninstall metapepview:
```Bash
pip uninstall metapepview
```