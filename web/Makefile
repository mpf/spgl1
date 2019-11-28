DOCS=$(patsubst %.jemdoc,%, $(wildcard *.jemdoc))
HDOCS=$(addsuffix .html, $(DOCS))

all : $(HDOCS)

%.html : %.jemdoc MENU spgl1.config
	jemdoc -c spgl1.config -o $@ $<

clean :
	-rm -f *.html eqs/*

update :
	@echo -n 'Copying to server...'
	rsync -av --delete-excluded \
		--exclude="*~"         \
		--exclude="Makefile"   \
		--exclude="*.config"   \
		--exclude="/*.jemdoc"  \
		--exclude=".DS_Store"  \
		--exclude="MENU"  \
                . \
                people.cs.ubc.ca:~mpf/public_html/spgl1
	@echo ' done.'
