set "comp=D:\Apps\MinGW\mingw64\bin\gfortran.exe"
set "objs=epcode_dep_qp.o epcode_dep_dp.o epcode_sub.o epcode_mod.o"
%comp% -O3 -c src\epcode_dep_qp.f
%comp% -O3 -c src\epcode_dep_dp.f
%comp% -O3 -c src\epcode_sub.f
%comp% -O3 -c src\epcode_mod.f
%comp% -O3 src\epcode_main_dp.f90 %objs% -o ep_dp.exe
%comp% -O3 src\epcode_main_qp.f90 %objs% -o ep_qp.exe
