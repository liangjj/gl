NAME=main
OBJECTS = $(NAME).f90
cc = ifort
MKLPATH=/opt/intel/composer_xe_2013.5.192/mkl/lib/intel64
MKLINCLUDE=/opt/intel/composer_xe_2013.5.192/mkl/include
MODPATH=/opt/intel/composer_xe_2013.5.192/mkl/include/intel64/lp64
FLAGS= -module $(MODPATH) \
 -assume byterecl -L$(MKLPATH) -I$(MKLINCLUDE)  -I$(MODPATH) -lmkl_lapack95_lp64 -lmkl_blas95_lp64  -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_intel_thread.a $(MKLPATH)/libmkl_core.a -Wl,--end-group -liomp5 -lpthread 

OUTNAME=$(NAME).out

$(NAME): $(OBJECTS)
	$(cc) -o $(OUTNAME) $(OBJECTS) $(FLAGS)
clean:
	rm -f $(OUTNAME) 

