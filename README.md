# g249 analysis

Set of codes to analyze the data of the g249 experiment.

## How to?

There are two types of programs to analyze the code `.cpp` and `.C` files. The first one have to be compiled using:

```bash
g++ setMultData.cpp $(root-config --cflags --libs)     -I$R3BROOT/r3bdata/footData     -I$R3BROOT/r3bdata/losData     -L$R3BROOT/../build/lib     -lR3BSsd     -lR3BLos     -o setMultData
```

Where `$R3BROOT` is the path where your R3BRoot installation is, and the just run the executable generated. On the other hand, `.C` files are ROOT macros that have to be run as:

```bash
root -l setMultData.C
``` 
