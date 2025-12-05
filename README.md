# ComTarget
ComTarget is a free, cross-platform, scientific software for  Chemical Target Prediction with Combinatorial Modeling. It is a command line oriented program. It can be executed in a shell (UNIX/Linux systems) or at a Windows or DOS command prompt. 
# Software Requirement
ComTarget requires the assistance of Obabel (Open Babel version 3.0.0 or higher) for format conversion and conformation search operations. Before using ComTarget, obabel must be installed in advance with environment variables properly configured, and it should be ensured that the confab module of obabel can be called normally. ComTarget also requires the assistance of AutoDock Vina to perform reverse docking operations. Please ensure that Vina has been installed and is callable.
## **Installation**
The ComTarget program is written in C and has been tested on Linux and Windows platforms. It is provided as source codes. To install the command line version of ComTarget please follow the instructions below.
1. Install GNU C compiler if you haven't already. On Ubuntu you may use sudo apt-get install gcc, for other Linux users you may need to use a different method.
2. Downloading the source code and info_data
3. compiling source code

```
gcc ComTarget.c -o ComTarget -lm -O3         
```
  or

```
gcc ComTarget.c -o ComTarget -lm 
```
If you are a Windows user, you can compile using gcc on the command prompt. Alternatively, install the free, portable, fast and simple C/C++ IDE Dev-C++, then Using Dev-C++, open the source code, click the compile option to compile and obtain the binary executable file.


## **How to Use**
Add ComTarget as an environment variable.
Open Terminal (Mac & Linux) or PowerShell (Windows).
The general synopsis for using ComTarget is:
```
ComTarget inputFile 
```

