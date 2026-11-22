# Installation

Meta-PepView is a python-based dashboard built with [Dash](https://dash.plotly.com/). It runs locally on your PC and can be accessed through the web browser. Before installation, check with the [[index#PC requirements|PC requirements]] if the system is suitable for running the dashboard tool. The recommended system requirements are:

- Operating system: meta-PepView was tested in Windows 11, Linux Ubuntu and Fedora. MacOS may work but is not actively supported.
- Minimal 16 GB RAM, 32 GB recommended for larger datasets.
- Minimal 10 GB free storage, additional space may be needed when working with large datasets.
- Web browser: Chromium-based browser, like Chrome or Edge, or Firefox.
- Internet connection required for downloading annotation databases and for Unipept annotation.

## Before installation: Software prerequisites

There are two ways to install meta-PepView: stand-alone installation using the pip package manager, or by creating a Docker container. 

Installing with `pip` requires Python to be installed. Meta-PepView can be installed on Python versions 3.11-3.15. Go to the [python website](https://www.python.org/downloads/) to download the installer, or for Linux, install python from the internal package manager. Check that a supported python version is installed (3.11-3.15). On Windows, make sure to enable "Add Python to PATH" during installation. To verify that Python is installed, open a terminal (or Command Prompt on Windows) and run:

```Bash
python --version
```

In addition to Python, Git has to be installed:

Windows: 
install the Git client [from their website](https://git-scm.com/install/windows). Note: Git initially only runs in a separate terminal applications, called 'Git Bash'. If Git is not active in Powershell or CMD, you may have to execute the installation instructions in Git Bash.

Ubuntu:
Git is usually already installed, but it may be manually Installed from their command line:

```Bash
sudo apt install git
```

Fedora:
Git is usually already installed, but it may be manually Installed from their command line:

```Bash
sudo dnf install git
```

For any operating system, to verify that Git is installed and active, run in the terminal: 
``` Bash
git --version
```

If meta-PepView is installed using a Docker container, [Docker](https://www.docker.com/) has to be installed. To run Docker containers in Windows, [WSL2](https://learn.microsoft.com/en-us/windows/wsl/) should be available, but this is usually already set up.

## Install meta-PepView with `pip`

The full meta-PepView application (including the [[build-reference-dataset|*mpv-buildref*]] command line utility) can be installed on any PC system with `pip`. This installation provides the meta-PepView dashboard, as well as a command line utility for construction of [[build-reference-dataset|experimental reference datasets]].

There are several ways to install meta-PepView this way: Directly with `pip`, or using package managers such as [pipx](https://pipx.pypa.io/stable/) or [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html). Installing with `pipx` or `conda` is recommended as it provides better isolation of meta-PepView from potential other Python. Below are the instructions to install meta-PepView using each method.

#### Installation with `pip`

When python (3.11-3.15) is installed, meta-PepView may be installed directly with `pip`. Installation via `pip` will add meta-PepView as package to the global python environment.

!!! note
    The instructions below execute `pip` directly from the command line. However if `pip` is not available, or if you want to run pip from a specific Python version, you may run `python -m pip ...`. Depending on installation method, you may have to change `python` to `py` or `python3`. For example: If you installed python using Microsoft Store, replace `python` with `python3`.

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

### Installation with `pipx`

First, make sure [python](https://www.python.org/) (3.11-3.15) is installed.
Once python is installed, install [pipx](https://pipx.pypa.io/stable/installation/) on the command line/terminal.

!!! note
    If multiple python interpreters are installed on the system, you can specify the interpreter to install the program on using the `--python` option: `pipx install --python 3.13 ...`. If you installed python using Microsoft Store, replace `py` with `python3` in the next line.

In Windows 10/11, `pipx` can be installed via `pip`:

```Powershell
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

### Installation inside a conda environment

The conda environment manager is an excellent tool to manage and isolate python programs. In case you want to manage projects using conda, it is recommended to install [miniconda](https://www.anaconda.com/docs/getting-started/concepts/anaconda-or-miniconda). This distribution provides a more minimal toolchain compared to the large Anaconda program. After installation, conda is often not yet exposed to the command line. The installer will inform you of this and provide steps to take. Usually, the following command should be executed:

```Bash
$ conda init
```

Close and reopen your terminal to activate the changes. You will see the current conda environment (usually "base") displayed on the command line. Once installed and setup, you may create new a environment using the following command:

```Bash
$ conda create -n environment_name python=3.15
```

Note that it is important to specify `python=...` in the command. Without it, conda will not install Python inside the environment. To ensure that the environment is suitable for installation of meta-PepView, the Python version should be between 3.11-3.15. Once the environment is created, it can be accessed by activating it:

```Bash
$ conda activate environment_name
```

The environment name should now be displayed on the command line. From here, you can simply install meta-PepView [[installation#Installation with `pip`|with pip]]. Note that meta-PepView tools can only be run when the conda environment is activated.

## Run meta-PepView inside a Docker container

To run the dashboard inside a container, make sure Docker is installed on the system.

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

Here, we start a container from the image "metapepview". The options`-p` publishes port 8050 from the container to port 8050 in the host pc, which allows us to reach the dashboard that is running inside the container. The dashboard can be accessed in the web browser at  `http://localhost:8050`. Note that the dashboard is only accessible from the host PC and not from other PC's in a local network.

## Download metaproteomics example data

Meta-PepView provides example datasets to help users explore the dashboard and test its features. However, these are not installed with meta-PepView, they have to be manually downloaded.

The example datasets are available in the [meta-PepView GitHub repository](https://github.com/RamonZwaan/metapepview/). Open the Meta-PepView GitHub repository, tutorial data is available in the Releases section on the right sidebar, under the name "tutorial datasets". You may download the data from there.

An instruction video that loads the metaproteomics datasets into the dashboard is present in [[getting-started|getting started]].
More details regarding data preparation and import are provided here: [[prepare-input-data|data prep]], [[data-import|data import]].
