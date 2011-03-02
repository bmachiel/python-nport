CONVERT=convert

EPSFILES := $(wildcard *.eps)
PNGFILES := $(patsubst %.eps,%.png,$(EPSFILES))

all: $(PNGFILES)

clean:
	-rm -rf $(PNGFILES)

%.png : %.eps
	$(CONVERT) $< $@
