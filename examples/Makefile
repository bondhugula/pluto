
DIRS = \
       covcol \
       corcol \
       dsyrk \
       dsyr2k \
       fdtd-2d \
       gemver \
       jacobi-1d-imper \
       jacobi-2d-imper \
       lu \
       matmul \
       mvt \
       seidel \
       ssymm \
       strmm \
       tmm

all: orig orig_par tiled par

orig:
	@-for d in $(DIRS); do \
		make -C $$d $@; \
		done

par:
	@-for d in $(DIRS); do \
		make -C $$d par; \
		done

pipepar:
	@-for d in $(DIRS); do \
		make -C $$d pipepar; \
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

