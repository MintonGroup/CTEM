### Installing CTEM for code development
************************************

CTEM uses CMake to automate the build process. You will need at least version 3.6.0.

Navigate to the topmost directory in your CTEM repository. It is best practice to create a ```build``` directory in your topmost directory from which you will compile CTEM. This way, temporary CMake files will not clutter up the project sub-directories. The commands to build the source code into a ```build``` directory and compile CTEM are:

```
$ cmake -B build -S .
$ cmake --build build
```
The CTEM executable, called `CTEM`, should now be created in the ```bin/``` directory.

If you wish to install CTEM into the system, execute the following command:

```
$ cmake --install build
```

This will install the project into the directory specified by `CMAKE_INSTALL_PREFIX` (the default is `/usr/local`), provided you have
sufficient permissions.


The CTEM CMake configuration comes with several customization options:

| Option                          | CMake command                                              |
| --------------------------------|------------------------------------------------------------|
| Build type                      | \-DCMAKE_BUILD_TYPE=[**RELEASE**\|DEBUG\|TESTING\|PROFILE] |
| Enable/Disable OpenMP support   | \-DUSE_OPENMP=[**ON**\|OFF]                                |
| Enable/Disable SIMD directives  | \-DUSE_SIMD=[**ON**\|OFF]                                  |
| Set Fortran compiler path       | \-DCMAKE_Fortran_COMPILER=/path/to/fortran/compiler        |
| Set path to make program        | \-DCMAKE_MAKE_PROGRAM=/path/to/make                        |
| Enable/Disable shared libraries (Intel only) | \-DBUILD_SHARED_LIBS=[**ON\|OFF]              |
| Add additional include path     | \-DCMAKE_Fortran_FLAGS="-I/path/to/libraries                 |
| Install prefix                  | \-DCMAKE_INSTALL_PREFIX=["/path/to/install"\|**"/usr/local"**] |


To see a list of all possible options available to CMake:
```
$ cmake -B build -S . -LA
```
For instance, if you wish to compile the project with debugging symbols enabled, run:

```
$ cmake -B build -S . -DCMAKE_BUILD_TYPE=Debug
```



The [CMake Fortran template](https://github.com/SethMMorton/cmake_fortran_template) comes with a script that can be used to clean out any build artifacts and start from scratch:

```
$ cmake -P distclean.cmake
```
