When contributing to open source projects, it's crucial to adhere to the developer documentation[1] to ensure compliance with the established coding standards and workflows. This alignment helps maintain the project's quality and facilitates effective collaboration among contributors.

Besides, we need to install Tamaas from source[2] to test our new branch.

#### C++ format

To obey the format to contribute this project, we can use extension wisely:

1. **Ensure `clang-format` is installed**:
   - On macOS or Linux, check if `clang-format` is installed by running `clang-format --version` in the terminal. If not installed, you can install it via Homebrew (on macOS) or a package manager (on Linux):
     ```bash
     brew install clang-format  # macOS
     sudo apt-get install clang-format  # Debian/Ubuntu Linux
     ```

2. **Configure VS Code to locate the `clang-format` executable**:
   - Open VS Code and go to Settings (`Cmd + ,` or `Ctrl + ,`).
   - Search for `clang-format.executable`.

Files can be formatted on-demand by right clicking in the document and selecting "Format Document", or by using the associated keyboard shortcut (usually Ctrl+⇧+F on Windows, Ctrl+⇧+I on Linux, and ⇧+⌥+F on macOS).

To automatically format a file on save, add the following to your vscode settings.json file:
```bash
{
    "editor.formatOnSave": true
}
```

###### *As redundancy measure, linting also occurs in continuous integration on Gitlab. Review the step artifacts to see which changes are necessary.*


## Installation from source, dependencies required for Tamaas:

- a C++ compiler with full C++14 and OpenMP support
- SCons (python build system)
- FFTW3
- thrust (1.9.2+)
- boost (pre-processor)
- python 3+ with numpy
- pybind11 (included as submodule)
- expolit (included as submodule)

Optional dependencies[2] are:

- an MPI implementation
- FFTW3 with MPI/threads/OpenMP (your pick) support
- scipy (for nonlinear solvers)
- uvw, h5py, netCDF4 (for dumpers)
- googletest and pytest (for tests)
- Doxygen and Sphinx (for documentation)

### Main steps for installation from source:

1. Clone from GitLab: 
```bash
git clone --recursive https://gitlab.com/tamaas/tamaas.git
```
2. Go to tamaas folder, and compile Tamaas with the default options by SCons
```bash
cd tamaas/
scons
```
#### SCons

Attention: To install SCons, which is a software construction tool used as a build system:

- Via **pip** (recommended for most users): Run `pip install scons`.
- Alternatively, if you prefer using a package manager on Linux: Use `sudo apt-get install scons` for Debian/Ubuntu or `sudo yum install scons` for Fedora/RHEL.

P.S. : If you don't have root access on the Linux machine to install SCons using `sudo apt-get install scons`, you can install SCons locally using Python's package manager, pip, which doesn't require root permissions:

- **Install SCons Locally Using Pip**:
     ```bash
     pip install --user scons
     ```
   This command installs SCons in your user directory, typically under `~/.local/bin`, which avoids the need for root privileges.

- **Update Your PATH**:
   - To use SCons from the terminal, ensure the directory where SCons was installed (`~/.local/bin`) is in your PATH. You can add this to your PATH by running:
     ```bash
     echo 'export PATH=$PATH:~/.local/bin' >> ~/.bashrc
     source ~/.bashrc
     ```

3. After compiling a first time, you can edit the compilation options in the file *build-setup.conf*


#### FFTW3

If you lack root access on a Linux machine, you can specify the file path in *build-setup.conf*.

#### Boost

If you lack root access on a Linux machine and can't install Boost via standard package managers, you can try a user-level installation or other methods. Here are some steps:

- **Install Boost from source**:
   - Download the Boost source code:
     ```bash
     wget https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.gz
     ```
   - Unpack the downloaded file:
     ```bash
     tar xfz boost_1_76_0.tar.gz
     ```
   - Go to the unpacked directory and configure with Bootstrap:
     ```bash
     cd boost_1_76_0
     ./bootstrap.sh --prefix=$HOME/local
     ```
   - Build and install:
     ```bash
     ./b2 install
     ```

   This will install Boost in the `local` folder in your home directory. Make sure your project build system can find the libraries in this directory.

- **Set environment variables**:
   - Add the `CPLUS_INCLUDE_PATH` environment variable to your `.bashrc` or `.bash_profile` to help the compiler locate the Boost headers:
     ```bash
     echo 'export CPLUS_INCLUDE_PATH=$HOME/local/include:$CPLUS_INCLUDE_PATH' >> ~/.bashrc
     source ~/.bashrc
     ```

#### thrust 
We can see more info by:
```bash
cat config.log 
```

For **thrust** and **cuda**, we can 

```bash
git clone https://github.com/NVIDIA/cccl.git
```

Automatically change *build-setup.conf* by:

```bash
scons THURST_ROOT=../cccl/thrust
scons LIBCUDACXX_ROOT=../cccl/libcudacxx
scons CXXFLAGS='-I/gagarine/temporaires/zli/cccl/libcudacxx/include
```

```bash
scons verbose=true
```

#### pybind11

```bash
git clone https://github.com/pybind/pybind11.git
scons PYBIND11_ROOT=../pybind11
```


```bash
scons build_tests=true
```

```bash
pip install -e build-release/pt
```



### We can also run tests with a virtual environment on your lab computer

First we merge master into our branch to get these changes.

You'll need to download a script to install pip

```bash
wget https://bootstrap.pypa.io/get-pip.py
```

Then create a virtualenv and boostrap pip

```bash
python3 -m venv --without-pip tamaas_venv
source tamaas_venv/bin/activate
python3 get-pip.py
pip install scons pytest
```

Then you can compile tamaas as usual, then run

```bash
pip install -e build-release/python
scons test
```

This should run a minimal number of tests. You can also `pip install scipy` if you wish, but we won't need scipy for the project.





#### Reference:

[1] https://tamaas.readthedocs.io/en/latest/developer.html

[2] https://tamaas.readthedocs.io/en/latest/quickstart.html#installation-from-source