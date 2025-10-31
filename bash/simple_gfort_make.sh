gfortran -O3 -o 1.o -c src/epcode_sub.f
gfortran -O3 -o 2.o -c src/epcode_dep_dp.f
gfortran -O3 -o 3.o -c src/epcode_dep_qp.f
gfortran -O3 -o 4.o -c src/epcode_mod.f
gfortran -O3 1.o 2.o 3.o 4.o -o ep_dp.exe src/epcode_main_dp.f90
gfortran -O3 1.o 2.o 3.o 4.o -o ep_qp.exe src/epcode_main_qp.f90
/bin/rm -f 1.o 2.o 3.o 4.o epcode_mod.mod 
