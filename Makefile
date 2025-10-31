########################################################################
#                                                                      #
#                       Makefile of the EP code                        #
#                                                                      #
########################################################################

FC      = gfortran
#FC      = ifort   

SRC     = src/
TST     = test/

FFLAGS  = -O3 -I$(SRC) 

DEPS = $(SRC)epcode_inc_dp.h  $(SRC)epcode_inc_qp.h

# 
# EP_main = \
#    ep_dp.f90 \
#    ep_qp.f90
# 
# EP_pack = \
#    epcode_sub.f      \
#    epcode_dep_dp.f   \
#    epcode_dep_qp.f   \
#    epcode_mod.f
#

EP_dp   = ep_dp.exe
EP_qp   = ep_qp.exe
EP_EXE  = $(EP_dp) $(EP_qp) 

EP_OBJ  = \
   epcode_sub.o      \
   epcode_dep_dp.o   \
   epcode_dep_qp.o   \
   epcode_mod.o

EP_OBJ_01 = \
   epcode_sub.o      \
   epcode_mod.o
EP_OBJ_dp = \
   epcode_dep_dp.o   
EP_OBJ_qp = \
   epcode_dep_qp.o 


EP_TEST = \
   ep_conf0_stdout.txt \
   ep_conf1_stdout.txt \
   ep_conf2_stdout.txt \
   ep_conf3_stdout.txt \
   ep_conf4_stdout.txt \
   ep_conf5_stdout.txt \
   ep_conf6_stdout.txt  

########################################################################

exe: $(EP_EXE) 
obj: $(EP_OBJ) 
all: $(EP_OBJ) $(EP_EXE) test  
test0: exe ep_conf0_stdout.txt 
test1: exe ep_conf1_stdout.txt 
test2: exe ep_conf2_stdout.txt 
test3: exe ep_conf3_stdout.txt 
test4: exe ep_conf4_stdout.txt 
test5: exe ep_conf5_stdout.txt 
test6: exe ep_conf6_stdout.txt 

%.o: $(SRC)%.f $(DEPS) 
	$(FC) $(FFLAGS) -o $@ -c $<

%_stdout.txt: $(TST)%.txt  
	./$(EP_dp) $^ > test_EPdp_$@ 
	./$(EP_qp) $^ > test_EPqp_$@ 

clean:
	@rm -rf $(EP_OBJ) epcode_mod.mod  test*.txt stdout*.txt

distclean:
	@rm -rf test*.txt stdout*.txt *.o *.mod *.exe  
 
#test: $(EP_dp) $(EP_qp) $(EP_TEST)
test: $(EP_TEST)
	@echo "All tests are done."

runEPdp: $(EP_dp) 
	./$(EP_dp) ep_conf.txt 

runEPqp: $(EP_qp) 
	./$(EP_qp) ep_conf.txt 

ep_dp.exe: $(EP_OBJ) $(SRC)/epcode_main_dp.f90 
	$(FC) $(FFLAGS) $(EP_OBJ) -o $(EP_dp) $(SRC)/epcode_main_dp.f90

ep_qp.exe: $(EP_OBJ) $(SRC)/epcode_main_qp.f90
	$(FC) $(FFLAGS) $(EP_OBJ) -o $(EP_qp) $(SRC)/epcode_main_qp.f90

########################################################################
