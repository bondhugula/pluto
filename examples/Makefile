
DIRS = ssymm \
	   strmm \
	   gemver \
	   tmm \
	   matmul \
	   dsyrk \
	   dsyr2k \
	   covcol \
	   corcol \
	   seidel \
	   jacobi-1d-imper \
	   jacobi-2d-imper \
	   fdtd-2d \
	   mvt \
	   lu \
	   heat-1dp \
	   heat-2dp \
	   heat-3dp

all: orig orig_par tiled par

orig:
	@-for d in $(DIRS); do \
		make -C $$d $@; \
		done

par:
	@-for d in $(DIRS); do \
		make -C $$d par; \
		done

lbpar:
	@-for d in $(DIRS); do \
		make -C $$d lbpar; \
		done

orig_par:
	@-for d in $(DIRS); do \
		make -C $$d $@; \
		done



tiled:
	@-for d in $(DIRS); do \
		make -C $$d tiled; \
		done

test:
	@-for d in $(DIRS); do \
		make -C $$d test; \
		done

perf: 
	@-for d in $(DIRS); do \
		echo "$$d"; \
		make --no-print-directory -C $$d perf; \
		done

clean:
	@-for d in $(DIRS); do \
		make -C $$d clean; \
		done

