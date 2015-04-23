
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
	   lu

all: orig orig_par tiled par distopt

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

dyn:
	@-for d in $(DIRS); do \
		make -C $$d dyn; \
		done

orig_par:
	@-for d in $(DIRS); do \
		make -C $$d $@; \
		done


tiled:
	@-for d in $(DIRS); do \
		make -C $$d tiled; \
		done

distomp:
	@-for d in $(DIRS); do \
		make -C $$d distomp; \
		done

distopt:
	@-for d in $(DIRS); do \
		make -C $$d distopt; \
		done

distopt_foifi: 
	@-for d in $(DIRS); do \
		make -C $$d distopt_foifi; \
		done

distopt_fop: 
	@-for d in $(DIRS); do \
		make -C $$d distopt_fop; \
		done

test:
	@-for d in $(DIRS); do \
		make -C $$d test; \
		done

dyn_test:
	@-for d in $(DIRS); do \
		make -C $$d $@; \
		done

dist_test: 
	@-for d in $(DIRS); do \
		make -C $$d dist_test; \
		done

dist_foifi_test: 
	@-for d in $(DIRS); do \
		make -C $$d dist_foifi_test; \
		done

dist_fop_test: 
	@-for d in $(DIRS); do \
		make -C $$d dist_fop_test; \
		done

dist_compare_test: 
	@-for d in $(DIRS); do \
		make -C $$d dist_compare_test; \
		done

dist_test2: 
	@-for d in $(DIRS); do \
		make -C $$d dist_test2; \
		done

dist_foifi_test2: 
	@-for d in $(DIRS); do \
		make -C $$d dist_foifi_test2; \
		done

dist_fop_test2: 
	@-for d in $(DIRS); do \
		make -C $$d dist_fop_test2; \
		done

dist_compare_test2: 
	@-for d in $(DIRS); do \
		make -C $$d dist_compare_test2; \
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

