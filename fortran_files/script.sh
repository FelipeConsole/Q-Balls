 #! /bin/csh

  # This is a simple shell script that will increment one step 
  # in the parameter that you choose to vary and passes it to the
  #  parameters.dat file that will be read by the program that calls colsys.f
  # to solve the system of BVP.

  # Be sure to cd into the right dir (not sure if needed)
  cd Desktop/Quali/Q-balls/fortran_files

  # compile the fortran code then execute it with the ./ option

for cf in $(seq 1.06 -0.01 0.00)
do
echo ${cf/,/.} > parametros.dat

cp solu_init_q_ball_ipr_415_xl_1200.dat solu.new

gfortran bolasq.f90 colsys.f -std=legacy -w
# ifort bolasq.f90 colsys.f -w
 
  ./a.out
  
  done
