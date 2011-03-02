CONVERT=convert
EPSTOPDF=epstopdf

EPSFILES := $(wildcard *.eps)
PNGFILES := $(patsubst %.eps,%.png,$(EPSFILES))
PDFFILES := $(patsubst %.eps,%.pdf,$(EPSFILES))

.PHONY: all clean png pdf

all: png pdf

clean:
	-rm -rf $(PNGFILES) $(PDFFILES)

png: $(PNGFILES)

pdf: $(PDFFILES)

%.png : %.eps
	$(CONVERT) $< $@

%.pdf : %.eps
	$(EPSTOPDF) $<
