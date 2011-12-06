##############################
# Gtest build.
##############################
# Build rules for libraries.
check_LTLIBRARIES = libgtest.la libgtest_main.la

libgtest_la_SOURCES = unittests/gtest/src/gtest-all.cc
libgtest_main_la_SOURCES = unittests/gtest/src/gtest_main.cc
libgtest_main_la_LIBADD = libgtest.la
libgtest_la_CXXFLAGS = -I$(top_srcdir)/unittests/gtest	\
-I$(top_srcdir)/unittests/gtest/include
libgtest_main_la_CXXFLAGS = -I$(top_srcdir)/unittests/gtest	\
-I$(top_srcdir)/unittests/gtest/include

# gtest source files that we don't compile directly.  They are
# #included by gtest-all.cc.
GTEST_SRC = unittests/gtest/src/gtest-death-test.cc	\
	    unittests/gtest/src/gtest-filepath.cc	\
	    unittests/gtest/src/gtest-internal-inl.h	\
	    unittests/gtest/src/gtest-port.cc		\
	    unittests/gtest/src/gtest-printers.cc	\
	    unittests/gtest/src/gtest-test-part.cc	\
	    unittests/gtest/src/gtest-typed-test.cc	\
	    unittests/gtest/src/gtest.cc

EXTRA_DIST = $(GTEST_SRC)

# Headers, not installed
GTEST_I = unittests/gtest/include/gtest
noinst_HEADERS += $(GTEST_I)/gtest-death-test.h				\
                  $(GTEST_I)/gtest-message.h				\
                  $(GTEST_I)/gtest-param-test.h				\
                  $(GTEST_I)/gtest-printers.h $(GTEST_I)/gtest-spi.h	\
                  $(GTEST_I)/gtest-test-part.h				\
                  $(GTEST_I)/gtest-typed-test.h $(GTEST_I)/gtest.h	\
                  $(GTEST_I)/gtest_pred_impl.h $(GTEST_I)/gtest_prod.h	\
                  $(GTEST_I)/internal/gtest-death-test-internal.h	\
                  $(GTEST_I)/internal/gtest-filepath.h			\
                  $(GTEST_I)/internal/gtest-internal.h			\
                  $(GTEST_I)/internal/gtest-linked_ptr.h		\
                  $(GTEST_I)/internal/gtest-param-util-generated.h	\
                  $(GTEST_I)/internal/gtest-param-util.h		\
                  $(GTEST_I)/internal/gtest-port.h			\
                  $(GTEST_I)/internal/gtest-string.h			\
                  $(GTEST_I)/internal/gtest-tuple.h			\
                  $(GTEST_I)/internal/gtest-type-util.h


