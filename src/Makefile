.PHONY: all fioss hminus

all:    fioss hminus
fioss:  fioss-64
hminus: hminus-64

OPTS="-j16"

clean-all clean-hminus clean-fioss: TARGET=clean

clean:        clean-all
clean-all:    all
clean-hminus: hminus
clean-fioss:  fioss

fioss-64: 
	make -f Makefile_fioss_64 $(TARGET) $(OPTS)

hminus-64:
	make -f Makefile_hminus_64 $(TARGET) $(OPTS)

