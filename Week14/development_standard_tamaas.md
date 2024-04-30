


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

{
    "editor.formatOnSave": true
}

###### *As redundancy measure, linting also occurs in continuous integration on Gitlab. Review the step artifacts to see which changes are necessary.*



#### Reference:

[1] https://tamaas.readthedocs.io/en/latest/developer.html
