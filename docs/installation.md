# Installation

Meta-PepView runs on top of the [Dash](https://dash.plotly.com/) framework in python. It is installed as a web server that is run locally on the PC and accessed inside the web browser. Before installation, check with the [[index#PC requirements|PC requirements]] if the system is suitable for running the dashboard tool.

Meta-PepView can be installed from `pip`/`pipx` or by starting a [Docker](https://www.docker.com/) container. 

### Install meta-PepView with `pipx` (or `pip`)

The full meta-PepView application (including the [[build-reference-dataset|*mpv-buildref*]] command line utility) can be installed on any PC system with `pip`. This installation provides the meta-PepView dashboard, as well as a command line utility for construction of [[build-reference-dataset|experimental reference datasets]]. 

It is recommended to install meta-PepView with [pipx](https://pipx.pypa.io/stable/). It works similar to `pip`, but is designed for installation of applications as opposed to libraries. When installing a package, `pipx` provides a separate environment for the tool to run without cluttering existing python environments and ensures the executables are exposed to the command line. When using `pip`, it is recommended to install meta-PepView in a separate environment ([virtual environment](https://docs.python.org/3/library/venv.html) or [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html))

!!! note
    In the instructions below, meta-PepView is installed by pip(x) by fetching the code from the GitHub repository. For this, git needs to be [installed on the PC](https://git-scm.com/install/linux).

#### Installation with `pipx`

First, make sure [python](https://www.python.org/) (3.11 or higher) is installed.
Once python is installed, install [pipx](https://pipx.pypa.io/stable/installation/) on the command line/terminal.

In Windows 10/11, `pipx` can be installed via `pip`:

```Powershell
# If you installed python using Microsoft Store, replace `py` with `python3` in the next line.
py -m pip install --user pipx
python -m pipx ensurepath
```

In Linux, `pipx` is often provided by a distro's built-in package manager:

Ubuntu:
```Bash
sudo apt install pipx
pipx ensurepath
```

Fedora:
```Bash
sudo dnf install pipx
pipx ensurepath
```

You may need to restart the terminal before you can run `pipx`.

Once pipx is installed, meta-PepView can be installed:
```Bash
pipx install git+https://github.com/RamonZwaan/metapepview.git
```

Once installed, the meta-PepView server can be started by running the application in the command line:
```Bash
metapepview
```

The dashboard can be accessed in the web browser (URL: `http://localhost:8050`).

To uninstall meta-PepView:
```Bash
pipx uninstall metapepview
```

#### Installation with `pip`

When python (3.11 or higher) is installed, meta-PepView may be installed directly with `pip`. Installation via `pip`
will add meta-PepView as package to the global python environment and it is not certain
that the application executables are exposed on the command line. To provide some isolation, it 
is recommended to install meta-PepView in a [virtual environment](https://docs.python.org/3/library/venv.html)
or in a [conda environment](https://docs.conda.io/en/latest/).

To install meta-PepView with pip:
```Bash
pip install git+https://github.com/RamonZwaan/metapepview.git
```

If the python environment is exposed to the operating system PATH, the meta-PepView server can be started by running the application in the command line:
```Bash
metapepview
```

To uninstall metapepview:
```Bash
pip uninstall metapepview
```

### Run meta-PepView inside a Docker container

To run the dashboard inside a container, make sure Docker (or another OCI compliant manager) is installed on the system.

First, build a [Docker image](https://docs.docker.com/get-started/docker-concepts/building-images/build-tag-and-publish-an-image/) using the Dockerfile template from the meta-PepView GitHub:

!!! note
    Make sure the Docker daemon is running during the command execution (start Docker desktop if installed).

```Bash
$ docker build -t metapepview https://github.com/RamonZwaan/metapepview.git#main
```

Now a Docker image named "metapepview" is created on the host system. The image can be found in Docker desktop, under *Images*.

Once the image is successfully built, start a container by running the following command:

```Bash
$ docker run -p 8050:8050 metapepview
```

Here, we start a container from the image "metapepview". The options`-p` publishes port 8050 from the container to port 8050 in the host pc, which allows us to reach the dashboard that is running inside the container. The dashboard can be accessed in the web browser at  `http://localhost:8050`. Note that the dashboard is only accessible from the host PC as long as the port is closed in the firewall.

