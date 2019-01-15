all: staticlib

rebuild: clean all

clean:
	$(RM) *.o libtinyvec3d.a

staticlib: libtinyvec3d.a

libtinyvec3d.a: tinyvec3d.o
	$(AR) $(ARFLAGS)s $@ $^

tinyvec3d.o: tinyvec3d.cpp tinyvec3d.h

