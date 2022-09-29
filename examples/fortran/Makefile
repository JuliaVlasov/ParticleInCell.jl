SRCS =	compute_rho_m6.f90 interpolation_m6.f90 meshfields.f90 particles.f90 \
	poisson_2d.f90 poisson_3d.f90 ua_steps.f90 ua_type.f90 fft.f90 \
        output.f90 interpolation_cic.f90 compute_rho_cic.f90

OBJS =	$(SRCS:.f90=.o)

LIBS =	-L/usr/local/lib -lfftw3 -lhdf5_fortran -lhdf5

F90 = gfortran

F90FLAGS = -I/usr/local/include -fstack-arrays -O0 -g -fbacktrace \
           -pedantic -Wuninitialized \
           -fcheck=array-temps,bounds,do,pointer,recursion \
           -ffpe-trap=invalid,zero,overflow

F90FLAGS = -I/usr/local/include -O3

LDFLAGS = $(F90FLAGS)  

PROG = bupdate efd test_poisson_3d test_pic_3d uapic3d
all: $(PROG)

bupdate: bupdate.f90 $(OBJS)
	$(F90) $(LDFLAGS) $< -o $@ $(OBJS) $(LIBS)

efd: efd.f90 $(OBJS)
	$(F90) $(LDFLAGS) $< -o $@ $(OBJS) $(LIBS)

test_poisson_3d: test_poisson_3d.f90 $(OBJS)
	$(F90) $(LDFLAGS) $< -o $@ $(OBJS) $(LIBS)

test_pic_3d: test_pic_3d.f90 $(OBJS)
	$(F90) $(LDFLAGS) $< -o $@ $(OBJS) $(LIBS)

uapic3d: uapic3d.f90 $(OBJS)
	$(F90) $(LDFLAGS) $< -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod *.o

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.mod.o:
	$(F90) $(F90FLAGS) -c $*.f90

meshfields.o: meshfields.f90
particles.o: particles.f90 meshfields.o
compute_rho_m6.o: compute_rho_m6.f90 meshfields.o particles.o ua_type.o
interpolation_m6.o: interpolation_m6.f90 meshfields.o particles.o compute_rho_m6.o
poisson_2d.o: poisson_2d.f90 meshfields.o 
ua_steps.o: ua_steps.f90 compute_rho_m6.o meshfields.o particles.o ua_type.o
ua_type.o: ua_type.f90 
fft.o: fft.f90
poisson_3d.o: poisson_3d.f90 meshfields.o
output.o: output.f90 meshfields.o
interpolation_cic.o: interpolation_cic.f90 particles.o
compute_rho_cic.o: compute_rho_cic.f90 particles.o
