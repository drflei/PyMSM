FC      = gfortran
FFLAGS  = -fno-automatic -ffree-line-length-none
LD      = $(FC) 
#LDFLAGS = -static
RM      = rm -f

SRCS = msm.for putils.for format.for

OBJS = $(SRCS:.for=.o)

MAIN = msm

$(MAIN): $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS) ./lib/liboneradesp_linux_x86_64

.phony: clean
clean:
	$(RM) $(OBJS) $(MAIN)

%.o: %.for interface.h 
	$(FC) -c $(IFLAGS) $(FFLAGS) $(FFSTAT) $< -o $@

