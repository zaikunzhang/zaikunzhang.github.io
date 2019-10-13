HOME=/home/zaikunzhang
LDIRECTORY=$(HOME)/Work/PersonalInformation/homepage
WWW=$(LDIRECTORY)/zaikunzhang.github.io
#JEMDOC=$(HOME)/local/bin/jemdoc
JEMDOC=$(HOME)/local/bin/jemdoc_mathjax

all: index blog software quotation coauthor conference misc powell teaching tools errors css makefile jemdoc $(JEMDOC)
	cp makefile $(WWW)/src

index: TARGET=index
index: index.jemdoc index.config MENU
	$(JEMDOC) -c $(TARGET).config $(TARGET)
	cp $(TARGET).html $(WWW)
	cp $(TARGET).* $(WWW)/src
	cp MENU $(WWW)/src
	cp $(TARGET).html ../

blog: TARGET=blog
blog: blog.jemdoc blog.config MENU
	$(JEMDOC) -c $(TARGET).config $(TARGET)
	cp $(TARGET).html $(WWW)
	cp -r eqs $(WWW)
	cp $(TARGET).* $(WWW)/src
	cp MENU $(WWW)/src
	cp -r eqs ../
	cp $(TARGET).html ../

software: TARGET=software
software: software.jemdoc software.config MENU
	$(JEMDOC) -c $(TARGET).config $(TARGET)
	cp $(TARGET).html $(WWW)
	cp -r eqs $(WWW)
	cp $(TARGET).* $(WWW)/src
	cp MENU $(WWW)/src
	cp $(TARGET).html ../
	cp -r eqs ../

coauthor: TARGET=coauthor
coauthor: coauthor.jemdoc coauthor.config MENU
	$(JEMDOC) -c $(TARGET).config $(TARGET)
	cp $(TARGET).html $(WWW)
	cp $(TARGET).* $(WWW)/src
	cp MENU $(WWW)/src
	cp $(TARGET).html ../

francais: TARGET=francais
francais: francais.jemdoc francais.config MENU
	$(JEMDOC) -c $(TARGET).config $(TARGET)
	cp $(TARGET).html $(WWW)
	cp $(TARGET).* $(WWW)/src
	cp MENU $(WWW)/src
	cp $(TARGET).html ../

quotation: TARGET=quotation
quotation: quotation.jemdoc quotation.config MENU
	$(JEMDOC) -c $(TARGET).config $(TARGET)
	cp $(TARGET).html $(WWW)
	cp $(TARGET).* $(WWW)/src
	cp MENU $(WWW)/src
	cp $(TARGET).html ../

conference: TARGET=conference
conference: conference.jemdoc conference.config MENU
	$(JEMDOC) -c $(TARGET).config $(TARGET)
	cp $(TARGET).html $(WWW)
	cp $(TARGET).* $(WWW)/src
	cp MENU $(WWW)/src
	cp $(TARGET).html ../

powell: TARGET=powell
powell: powell.jemdoc powell.config MENU
	$(JEMDOC) -c $(TARGET).config $(TARGET)
	cp $(TARGET).html $(WWW)
	cp $(TARGET).* $(WWW)/src
	cp MENU $(WWW)/src
	cp $(TARGET).html ../

misc: TARGET=misc
misc: misc.jemdoc misc.config MENU
	$(JEMDOC) -c $(TARGET).config $(TARGET)
	cp $(TARGET).html $(WWW)
	cp $(TARGET).* $(WWW)/src
	cp MENU $(WWW)/src
	cp $(TARGET).html ../

teaching: TARGET=teaching
teaching: teaching.jemdoc teaching.config MENU
	$(JEMDOC) -c $(TARGET).config $(TARGET)
	cp $(TARGET).html $(WWW)
	cp $(TARGET).* $(WWW)/src
	cp MENU $(WWW)/src
	cp $(TARGET).html ../

tools: TARGET=tools
tools: tools.jemdoc tools.config MENU
	$(JEMDOC) -c $(TARGET).config $(TARGET)
	cp $(TARGET).html $(WWW)
	cp $(TARGET).* $(WWW)/src
	cp MENU $(WWW)/src
	cp $(TARGET).html ../

errors: errors/40*.jemdoc errors/errors.config 403 404 MENU

403: TARGET=403
403: 
	cd ./errors && $(JEMDOC) -c errors.config $(TARGET)
	if test -d $(WWW); \
	then echo $(WWW) exists.; \
	else echo "mounting remote directory ..."; \
	mkdir -p $(LDIRECTORY); sshfs $(RDIRECTORY) $(LDIRECTORY); \
	fi
	cp errors/$(TARGET).html $(WWW)/errors
	cp errors/$(TARGET).* $(WWW)/src/errors
	cp errors/*.config $(WWW)/src/errors
	cp MENU $(WWW)/src
	cp errors/$(TARGET).html ../errors
404: TARGET=404
404:
	cd ./errors && $(JEMDOC) -c errors.config $(TARGET)
	cp errors/$(TARGET).html $(WWW)/errors
	cp errors/$(TARGET).* $(WWW)/src/errors
	cp errors/*.config $(WWW)/src/errors
	cp MENU $(WWW)/src
	cp errors/$(TARGET).html ../errors

css: jemdoc.css
	cp jemdoc.css zhangzk.css
	cp *.css ../
	cp *.css $(WWW)
	cp *.css $(WWW)/src

jemdoc: $(JEMDOC)
	cp $(JEMDOC) ./
	cp $(JEMDOC) ../
	cp $(JEMDOC) $(WWW) 
	cp $(JEMDOC) $(WWW)/src

htaccess:
	cp .htaccess ../
	cp .htaccess $(WWW) 
	cp .htaccess $(WWW)/src

sitemap:
	# generate sitemap at https://www.xml-sitemaps.com
	cp sitemap* urllist.txt ror.xml ../
	cp sitemap* urllist.txt ror.xml $(WWW)
	cp sitemap* urllist.txt ror.xml $(WWW)/src

robots:
	cp robots.txt ../
	cp robots.txt $(WWW)
	cp robots.txt $(WWW)/src

clean: 
	rm -f index.html blog.html software.html quotation.html francais.html errors/40*.html
