if HAVE_TEXI2DVI
cloog.pdf: cloog.texi gitversion.texi
	$(TEXI2DVI) --pdf $<
endif

gitversion.texi: @GIT_INDEX@
	echo '@set VERSION '`$(top_builddir)/genversion.sh`'' > $@
