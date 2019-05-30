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

### Setup on the Cluster

The AerOpt GUI depends on a specific version of the AerOpt tool on the cluster.
You can use git to get the newest version.
Log on to the cluster and run
```
mv AerOpt AerOpt-old
git clone https://github.com/DrBenEvans/AerOpt.git
```

Compile the new version of AerOpt
```
module load mpi
cd AerOpt/source
make all
ls -s AerOpt_v3.5 ../AerOpt
```

You should now be able to run jobs on the cluster using the GUI.


## Deploying

Deploying is essentially creating a compiled executable you can distribute without needing any additional libraries or tools.
There are two deployed releases in [Jarno's github page](https://github.com/rantahar/AerOpt-gui/releases).
All you need to do to run this version of AerOpt is to donwload AerOptGui.zip, extract the file and run AerOptGui.exe in the new folder.
To be able to run on a cluster, you also need to download AerOpt.tar.gz
on the cluster.
```
mv AerOpt AerOpt_backup
wget https://github.com/rantahar/AerOpt-gui/releases/download/0.2/AerOpt.tar.gz
tar -xvf AerOpt.tar.gz
```

This is much easier for a user than compiling the whole thing.

To create a new release compile the program using `release` setting in QT Creator.
This is above the play symbol on the bottom left in QT Creator 4.9.
When you compile the program, an executable called `AerOptGui.exe` is created in the same directory the AerOptGui folder is in.

Download the [zip archive](https://github.com/rantahar/AerOpt-gui/releases/download/0.2/AerOptGui.zip) and extract it.
The folder contains everything AerOptGui needs to run.
Replace `AerOptGui.exe` with the newly compiled executable and it should work directly.
