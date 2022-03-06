# Q-Balls

Q-balls are an example of solitons.
Solitons ares ...

More spcefically, Q-balls are non-topological solitons. Obviously, this means that they are not topological.
Which means that there are topological solitons.

Topological solitons are ...


The Q-ablls I am describinh here was first introduced by Coleman in 1985 paper [...]


In this repository I included the files necessary to solve the Q-balls equations

This is some topic I wou;d like to disccuss further. 

These equatins, are non linear. Whci makes things extremelhy difficult.
In special, THEY ARE SINGULAR AT THE ORIGIN (a SINGULAR REGULAR POISN, NONTHELES).

tHE BEST SOLVER i COULD FIND IS THIS 1980 FORTRAN SOLVER CALLD COLSYS.

Colsys uses an adaptive grid selection and a newton rapshson method to solve the boundary value Problems.

I tried several options to solve the equations (Mathematica, Julia, Matlab...) but none of them worked as good as Colsys.
There is a draopback thoug, It's not as easy to use as Mathematica for example.


Now for the contents in this repo.

I included the colsys program (the BVP solver) and the program that calls colsys with approprtiate equtions to be solved.


There ś also a Mathemativa notebook included to derive the equations, and a bash script to automazi the execution of the program for different values of the input parameters.  

