.PHONY: all fioss hminus fioss-32 fioss-64 hminus-32 hminus-64
all:    fioss hminus
32:     fioss-32 hminus-32
64:     fioss-64 hminus-64
fioss:  fioss-32 fioss-64
hminus: hminus-32 hminus-64
#OPTS="-j4"
OPTS="-j16"
all-clean hminus-clean hminus-32-clean hminus-64-clean: TARGET=clean
32-clean 64-clean fioss-clean fioss-32-clean fioss-64-clean: TARGET=clean



clean-all clean-hminus clean-hminus-32 clean-hminus-64: TARGET=clean
clean-32 clean-64 clean-fioss  clean-fioss-32 clean-fioss-64: TARGET=clean

clean: clean-all
clean-all: all
clean-64: 64
clean-32: 32
clean-hminus: hminus
clean-hminus-32: hminus-32
clean-hminus-64: hminus-64
clean-fioss: fioss
clean-fioss-32: fioss-32
clean-fioss-64: fioss-64

fioss-32:
	make -f Makefile_fioss_32 $(TARGET) $(OPTS)

fioss-64: 
	make -f Makefile_fioss_64 $(TARGET) $(OPTS)

hminus-32:
	make -f Makefile_hminus_32 $(TARGET) $(OPTS)

hminus-64:
	make -f Makefile_hminus_64 $(TARGET) $(OPTS)

