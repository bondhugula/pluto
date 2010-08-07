noinst_PROGRAMS = @examples@
EXTRA_PROGRAMS = example32 example64 exampleMP

example32_SOURCES = example.c
example32_CPPFLAGS = @DFLAGS32@ -I$(top_srcdir)/include
example32_LDADD = ../libpiplib32.la
example64_SOURCES = example.c
example64_CPPFLAGS = @DFLAGS64@ -I$(top_srcdir)/include
example64_LDADD = ../libpiplib64.la
exampleMP_SOURCES = example.c
exampleMP_CPPFLAGS = @MP_DFLAGS@ -I$(top_srcdir)/include
exampleMP_LDADD = ../libpiplibMP.la @EXTRA_LIBS@

PIPTEST = \
	big \
	boulet \
	brisebarre \
	cg1 \
	esced \
	ex \
	ex2 \
	expansion \
	fimmel \
	max \
	negative \
	small \
	sor1d \
	square \
	sven

EXTRA_DIST = \
	$(PIPTEST:%=%.pip) \
	$(PIPTEST:%=%.ll)

test:
	@failedtest=0; \
	for x in $(PIPTEST) ; do \
		echo "Verify file $(srcdir)/$$x" ; \
		./$(EXAMPLE) < $(srcdir)/$$x.pip > xyz ;\
		diff -w xyz $(srcdir)/$$x.ll ; \
		result=$$?; \
		if [ "$$result" -ne "0" ]; then \
			echo "Error: $$x is not the same"; \
			failedtest=`expr $$failedtest + 1`; \
		else \
			echo "$$x passed"; \
	        fi; \
	done ; \
	rm xyz ; \
	if [ $$failedtest != 0 ]; then \
		echo "$$failedtest tests failed"; \
	else \
	        echo "PIP works correctly :-)"; \
        fi
