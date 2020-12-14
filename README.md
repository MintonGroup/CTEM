### Installing CTEM for code development
************************************

CTEM uses Autotools to automate makefile generation. To set up your system, 
first run the script 

`$ ./autogen.sh`

You may have to modify the contents of src/Makefile.am in order to set the 
appropriate flags for your compiler.  The default file contains examples of 
common flags (optimized code vs. debugging code) for gfortran and ifort.  Simply
activate the appropriate set of flags for your system.

From then on, you simply need to execute the following:

```
$ cd build
$ ../configure
$ make clean
$ make
```

The executable that is generated is `build/src/CTEM`


************************************
### Running CTEM
************************************
An example CTEM simulation is included in the idl directory. You can execute the
CTEM run using the IDL frontend using:

`$ idl < start.in`

Alternatively you may wish to use the Python front end.


`$ python3 ctem_driver.py`

The run parameters are controlled by the ctem.in file.
