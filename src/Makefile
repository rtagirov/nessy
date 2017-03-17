.PHONY: all fioss hminus

all:    FIOSS HMINUS
fioss:  FIOSS
hminus: HMINUS

OPTS="-j16"

clean-all clean-hminus clean-fioss: TARGET=clean

clean:        clean-all
clean-all:    all
clean-hminus: hminus
clean-fioss:  fioss

FIOSS: 
	make -f Makefile_fioss $(TARGET) $(OPTS)

HMINUS:
	make -f Makefile_hminus $(TARGET) $(OPTS)
