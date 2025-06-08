check_PROGRAMS = picopic_unit

picopic_unit_SOURCES = unit/grid3dTest.cpp\
						unit/constantTest.cpp\
						unit/geometryTest.cpp\
						unit/gridTest.cpp\
						unit/msgTest.cpp\
						unit/vector3dTest.cpp\
						unit/main.cpp

picopic_unit_LDFLAGS = -lgtest -lgmock

# math_interp_SOURCES = math/interp.cpp
# math_interp_LDFLAGS = -pthread -lgtest -lgmock

# math_flowing_timeseries_SOURCES = math/flowing_timeseries.cpp
# math_flowing_timeseries_LDFLAGS = -pthread -lgtest -lgmock

TESTS = $(check_PROGRAMS)

functional-test:
	@python3 functional/test.py

perf-test:
	@python3 performance/test.py -l 5

check-local: functional-test perf-test
