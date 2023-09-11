### Installing CTEM for code development
************************************

CTEM uses CMake to automate the build process. You will need at least version 3.6.0.

To build the executable and all artifacts in a directory called `build/`, simply run the following commands from the CTEM root directory:

```
$ cmake -B build -S .
$ cmake --build build
```

To clean all build artifacts, you can execute the included `distclean.cmake` script with the following command:

```
$ cmake -P distclean.cmake
```