DOCS=$(patsubst %.jemdoc,%, $(wildcard *.jemdoc))
HDOCS=$(addsuffix .html, $(DOCS))

all : $(HDOCS)

%.html : %.jemdoc MENU spgl1.config
	jemdoc -c spgl1.config -o $@ $<

clean :
	-rm -f *.html eqs/*

