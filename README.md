# The AerOpt GUI

## Compiling in QT Creator

### Install Compilation Requirements

The recommended way of compiling the AerOpt GUI is using the QT Creator IDE.
You can download the IDE on the [QT website](https://www.qt.io/download).
Choose the open source option and download the installer.
In the installer, you will need to create a QT account if you don't already have one.
When selecting components, check the latest version of QT and 64 bit version of MinGW.

You will also need the ssh and ssl libraries.
On Windows, we recommend using [vcpkg](https://github.com/microsoft/vcpkg) for this.
First you will need to install [git](https://gitforwindows.org/), if you don't already have it.
Go to the [git for Windows home page](https://gitforwindows.org/) and click download.

After installing git, you should also have a program called `Git Bash`.
Run it to open a command line interface. If you are familiar with it, use `cd`
to navigate where you want to install vcpkg.
If not, you probably want to install it in your home folder, which is where you start.

Copy the following commands in to the command line and press Enter to execute
```
git clone https://github.com/Microsoft/vcpkg.git
cd vcpkg
.\bootstrap-vcpkg.sh
```
For later convenience, you may want to execute the following line as well.
Windows will ask for admin access.
```
./vcpkg integrate install
```

Now you can install libraries. You are already in the correct folder, but once you
close the window and open a new one, you will need to run
```
cd vcpkg
```
to get back.
Now install `libssh`, `openssl` and `zlib` with the following command
```
./vcpkg install libssh:x64-windows openssl:x64-windows zlib:x64-windows
```

Next, open `System Properties` in the start menu, and press enter.
Click `Environment Variables`, find `path` in the second panel and double click.
Add the following three lines as separate entries
```
C:\Users\USERNAME\vcpkg\packages\libssh_x64-windows\bin
C:\Users\USERNAME\vcpkg\packages\openssl-windows_x64-windows\bin
C:\Users\USERNAME\vcpkg\packages\zlib_x64-windows\bin
```
replacing `USERNAME` with your actual user name. For this change to take effect, you may need to reboot.

Note that if you installed `vcpkg` in a custom location, you need to change the paths accordingly.

### Setup Project

Next, download the AerOpt GUI. Since you have git, you can use it to clone the project. Open `Git Bash` and run the command
```
git clone https://github.com/DrBenEvans/AerOptGUI.git
```

Open QT Creator, click `Open Project` and open `AerOptGUI\GeneticGui.pro`, which should now be in your home directory.
The first time you do this the IDE will ask you about the configuration.
Choose the latest 64 bit MinGW version and click `configure`.

Next, open `GeneticGui.pro` in the editor and scroll to the bottom.
You will see several lines including the path
```
C:\Users\USERNAME\src\vcpkg...
```
Replace `USERNAME` with your username on each line. If `vcpkg` is in a custom location, use that path instead.