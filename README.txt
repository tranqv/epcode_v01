We are now in the EPCODE directory. Suppose that our computer operating system is Linux. 

1. Performing the 'ls | nl' command we should see

      1  bash
      2  batch
      3  epconf_default.txt
      4  Makefile
      5  README
      6  src
      7  test

   The following entries are directories:

      bash, batch, src, test 

2. The directory 'src' contains all routines of EPCODE. Upon executing 
   the command 'ls src/ | nl' we should see 

      1  epcode_dep_dp.f
      2  epcode_dep_qp.f
      3  epcode_inc_dp.h
      4  epcode_inc_qp.h
      5  epcode_main_dp.f90
      6  epcode_main_qp.f90
      7  epcode_mod.f
      8  epcode_sub.f

   where 

      +  "_dp", "_qp":
            Suffixes indicating the double and quadruple precisions for real numbers, 
            respectively. 
      +  "epcode_dep_dp.f", "epcode_dep_qp.f":
            Files including subroutines/functions selected from ARPACK, LAPACK, and 
            BLAS, which were modified slightly to support specified precision for 
            numbers. These routines are required by EPCODE.
      +  "epcode_inc_dp.h" and "epcode_inc_qp.h" define working precisions.
      +  "epcode_main_dp.f90", "epcode_main_qp.f90":
            Main programs which read a configuration file to specfify the values
            of parameters and call the core of EPCODE to perform all calculations.
      +  "epcode_mod.f":
            Module defines generic interface for various EPCODE subroutines, which 
            serve for different scenarios of parameter input.
      +  "epcode_sub.f":
            The core of EPCODE. 

3. The directory 'test' includes seven examples of configuration file of EPCODE. 
   Performing the command 'ls test/ | nl' we should see 

      1   ep_conf0.txt   - [Format A]
      2   ep_conf1.txt   - [Format A]	
      3   ep_conf2.txt   - [Format A]
      4   ep_conf3.txt   - [Format A]
      5   ep_conf4.txt   - [Format A]
      6   ep_conf5.txt   - [Format B] 
      7   ep_conf6.txt   - [Format B] 


    Formats:
    - Format A: For deformed nuclei with non-degenerate single-particle energies.
    - Format B: For spherical nuclei with degenerate single-particle energies.

	Here, "ep_conf0.txt" provides guidance for writing in Format A for deformed nuclei,
	while "ep_conf5.txt" offers instructions for Format B used with spherical nuclei having degenerate energy levels.
	
	For more details:
	+  Descriptions of the declarations are provided directly within the configuration files.
	
	+  If we want to comment inside a configuration file, use "#" or "!".
   		--------------------------------
   		- Example: ep_conf1.txt
   		--------------------------------
   

	+  To check whether the current computer configuration meets the computational requirements for the given values of nome and npar,
   	   the user can enable the "prep = true" flag. 
   		--------------------------------
   		- Example: ep_conf2.txt
   		--------------------------------
   

	+  To define G, i.e., 'g', for the case that G is given as a matrix, check 
   		--------------------------------
  		- Example: ep_conf3.txt 
  		--------------------------------
   		- We define G(:,:) element by element. Each line provides one entry in the format: i    j    G(i,j).
   

	+  To define the single-particle spectrum E(:), i.e., 'e', we must provide all elements either in a single row or arranged vertically in a single column.
	   If you intend to define it across multiple lines, leave the right-hand side of the = sign empty. Otherwise, the program may not work properly.
   		--------------------------------
   		- Example: ep_conf4.txt
   		--------------------------------
   		- If we do not define 'e', the program set E(k) = k, for k = 1,NOME. 
   
   
	+  To perform calculations for a spherical nucleus, we provide a test case using the realistic nucleus 114Sn, 
   	   incorporating experimental single-particle energies and pairing matrix elements V0 (i    j    V_0(ii,jj)), check:
   		--------------------------------
   		- Example: ep_conf5.txt
   		--------------------------------
   		The code will automatically calculate G_{ij} from V_0(ii,jj) using the formula: G_{ij} = 2 V_0(ii,jj) / [(2i+1)(2j+1)]^(1/2)
   
 	+ To perform calculations for a spherical nucleus using a scalar G, check:
   		--------------------------------
   		- Example: ep_conf6.txt
   		--------------------------------

3. The directories "bash" and "batch" provide Linux and Windows users with 
   Bash and Batch scripts, respectively, to compile EPCODE.

4. "ep_conf.txt" is the default configuration file, which will be used 
   when we run EPCODE executable program alone without argument.

5. "Makefile" is Linux makefile which defines conventional commands:
   +  make 
         to compile the whole of EPCODE and run all the tests.
   +  make clean 
         to clean all output files, except executable EPCODE files.
   +  make distclean 
         to clean all output files.
   +  make test 
         to run all the tests after we compile EPCODE
   In the default setting, the makefile works with Gfotran compiler. 
   Users can edit FC in Makefile to work with other Fortran compilers.
   
6. If the user does not want to use Makefile but just compile and run manually, then perform the following steps:  

   +	Step 1: run these commands: 
   
           gfortran -O3 -o 1.o -c src/epcode_sub.f
           gfortran -O3 -o 2.o -c src/epcode_dep_dp.f
           gfortran -O3 -o 3.o -c src/epcode_dep_qp.f
           gfortran -O3 -o 4.o -c src/epcode_mod.f
           
   +	Step 2: run 
   
           gfortran -O3 1.o 2.o 3.o 4.o -o ep_dp.exe src/epcode_main_dp.f90
           
           to create the ep_dp.exe file with double precision
           or run 
           
           gfortran -O3 1.o 2.o 3.o 4.o -o ep_qp.exe src/epcode_main_qp.f90
           
           to create the ep_qp.exe file with quadruple precision

   +	Step 3: run exe file with config file 
   
           ./ep_dp.exe test/ep_conf1.txt
           
           or
           
           ./ep_qp.exe ep_conf.txt

 Good luck and enjoy.
