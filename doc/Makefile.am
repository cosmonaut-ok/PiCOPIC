#
# doxygen
#
if HAVE_DOXYGEN
directory = $(top_srcdir)/doc/man/man3

man_MANS = $(directory)

$(directory): doxyfile.stamp

doxyfile.stamp: Doxyfile
	mkdir -p img
	for i in $(shell find . -name '*.svg');\
	do if [ ! -f img/$$(dirname $$i)/$$(basename $$i .svg).png ];\
		then echo "converting $$i to img/$$(dirname $$i)/$$(basename $$i .svg).png";\
		convert $$i img/$$(dirname $$i)/$$(basename $$i .svg).png 2>/dev/null;\
		fi;\
	done
	$(DOXYGEN) $^
	$(MAKE) -C out/latex
	echo Timestamp > $@

CLEANFILES = doxyfile.stamp

all-local: doxyfile.stamp
clean-local:
	-rm -rf $(top_srcdir)/doc/out $(top_srcdir)/doc/doxyfile.stamp $(top_srcdir)/img

endif
