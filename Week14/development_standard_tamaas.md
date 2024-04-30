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

#### SCons

To install SCons, which is a software construction tool used as a build system:

- Via **pip** (recommended for most users): Run `pip install scons`.
- Alternatively, if you prefer using a package manager on Linux: Use `sudo apt-get install scons` for Debian/Ubuntu or `sudo yum install scons` for Fedora/RHEL.

#### thrust 




#### Reference:

[1] https://tamaas.readthedocs.io/en/latest/developer.html

[2] https://tamaas.readthedocs.io/en/latest/quickstart.html#installation-from-source