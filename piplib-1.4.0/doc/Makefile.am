EXTRA_DIST = images piplib.texi

if HAVE_TEXI2DVI
piplib.pdf: piplib.texi gitversion.texi
	$(TEXI2DVI) --pdf $<
endif

gitversion.texi: @GIT_INDEX@
	echo '@set VERSION '`$(top_builddir)/genversion.sh`'' > $@
