#FC = mpif90 -O3 -ipo -xHOST -openmp ##T2K-Tsukuba@Tsukuba w openMP
#FC = mpif90 -O3 -ipo -xHOST ##T2K-Tsukuba@Tsukuba w/o openMP
#FC = ifort -O3 -ipo -xHOST  ##kashiwa@ISSP w/o openMP
#FC = mpifrt -Kfast -Fixed -X9 ##BX900@JAEA w/o openMP (maybe)
#FC = mpif90 -O3 -ipo -xHOST -openmp ##cal*@nucl_thr w openMP
#FC = mpif90 -O3 -ipo -xHOST ##cal*@nucl_thr w/o openMP
#FC = mpif90 -O3 -mkl  ##cal*@nucl_thr w/o openMP
#FC = mpif90 -O3 -ip -ipo -xHOST  -mkl  ##cal*@nucl_thr w/o openMP
FC = mpif90 -O3 -ip -ipo -xHOST  -mkl ##cal*@nucl_thr w/o openMP
#FC = mpif90 -O0 -mkl ##cal*@nucl_thr w/o openMP
#FC = mpif90 -O3 -ipo -xHOST -mkl ##cal*@nucl_thr w/o openMP
#FC = mpif90 -O3 -ipo -xHOST  ##cal*@nucl_thr w/o openMP
#FC = mpifrtpx -Kfast,openmp ##K@AICS

#LN = -lmpi #kashiwa
LN = #other
#LN = -lmkl_lapack -lmkl_em64t #typical intel MKL including LAPACK, fftw
#LN = -mkl=parallel #cal4, cal7, gpu1, gpu2: modern intel MKL


VPATH = src:object
SRC = $(shell cd src ;ls *.f90 ;cd ..)
OBJ = $(SRC:.f90=.o)
OBJ_dir = $(addprefix object/,$(OBJ))

PROG = crab

$(PROG):global_variables.o PSE_variables.o main.o PSE_band_calc_variables.o $(OBJ)
	$(FC) -o $(PROG) $(OBJ_dir) $(LN)

main.o:main.f90
	$(FC) -c $< $(LN);mv $@  object 
%.o:%.f90
	$(FC) -c $< $(LN);mv $@  object 


clean:
	rm  -f  object/*.o  *.mod crab
clean_complete:
	rm  -f *~  */*~ */*/*~ object/*.o  */*.mod *.mod crab */#*
