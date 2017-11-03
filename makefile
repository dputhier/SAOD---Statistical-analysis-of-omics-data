################################################################
## Synchronize web site on the server

MAKEFILE=makefile

help:
	@echo 
	@echo "-Avalaible targets:"
	@perl -ne 'if (/^([a-z]\S+):/){ print "\t$$1\n";  }' ${MAKEFILE}
	@echo
	@echo "-Some examples:"
	@echo "\tmake checkRCodeInHtml HTMLFILE=/path/to/html/file.html"
	@echo "\tmake extractRCodeFromHtml HTMLFILE=/path/to/html/file.html"
	@echo

USER=teacher
HOST=pedagogix-tagc.univ-mrs.fr
WDIR=/homes/teacher/courses/ASG1



## Directory to synchronize
EXCLUDED=--exclude jobs  \
	--exclude '*~' \
	--exclude '*.old' \
	--exclude '*.log' \
	--exclude '.DS_Store' \
	--exclude '.RData' \
	--exclude '.RHistory' \
	--exclude 'screen_pictures_*' \
	--exclude '*.Rproj.user' \
	--exclude '*.RData' \
	--exclude '*.Rhistory' \
	--exclude annotation_projects
RSYNC=rsync -ruptvl -e ssh -z  ${EXCLUDED} 

################################################################
## Clean temporary files created by emacs
clean:
	find . -name '*~' -exec rm {} \;
	find . -name '.#*' -exec rm {} \;
	find . -name '.DS_Store' -exec rm {} \;

################################################################
## Publish on the web site

TO_SYNC=*
publish: clean
	@echo "Transfering as $(USER) to $(HOST) in wdir=$(WDIR)"
	@rsync -ruptvl ${EXCLUDED} ${OPT} ${TO_SYNC}  $(USER)@$(HOST):$(WDIR)


pushToGit:
	@git pull
	@git add .
	@git commit -m "Added something"
	@git push

pushAndPublish: pushToGit publish


################################################################
##  synchronize the slides from Jacques' folder stats4bioinfo
sync_slides:
	rsync -ruptvl  \
		--exclude '*6spp.pdf' \
		--exclude "old*" \
		--exclude "prerequisites.pdf" \
		--exclude "*exercise*" \
		--exclude "statistiques_bioinformatique.pdf" \
		~/statistics_bioinformatics/pdf_files/*.pdf pdf_files/

################################################################
## Browse the Web site
#OPEN=open -a ${BROWSER}
#BROWSER=firefox
OPEN=open
local:
	${OPEN} http://localhost/courses/ASG/

web_github:
	${OPEN} https://dputhier.github.io/ASG/

## The Web site was previously hosted on pedagogix. We should move it
## but there are still some files to be transferred.
web_tagc:
	${OPEN} http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/

################################################################
# Check R Code.
# Note that R code should be embeded in the following markup:
# <pre class="brush:r">  # or <pre class='brush:r'>
#  R code
# </pre>

ifdef MAKECMDGOALS
 ifeq ($(MAKECMDGOALS), checkRCodeInHtml)
  ifeq ($(HTMLFILE),)
   $(error HTMLFILE variable is undefined. Use: "make checkRCodeInHtml HTMLFILE=/path/to/html/file.html")<
  endif
 endif
endif

# This temporary file store the html code from which comments have been disgarded.
TMPFILE:=$(shell mktemp /tmp/check.XXXXXX)

checkRCodeInHtml:
	@perl -ne 'print unless (/<!--/../-->/)' $(HTMLFILE)  > $(TMPFILE)
	@perl -ne 'print if(/<pre\s+class=\s*[\047"]brush:r[\047"]\s*>/../<\/pre>/)' $(TMPFILE) | grep -Ev "<.?pre"  | R --vanilla
	@rm -f $(TMPFILE)

ifdef MAKECMDGOALS
 ifeq ($(MAKECMDGOALS), extractRCodeFromHtml)
  ifeq ($(HTMLFILE),)
   $(error HTMLFILE variable is undefined. Use: "make extractRCodeFromHtml HTMLFILE=/path/to/html/file.html")<
  endif
 endif
endif


extractRCodeFromHtml:
	@perl -ne 'print unless (/<!--/../-->/)' $(HTMLFILE)  > $(TMPFILE)
	@perl -ne 'print if(/<pre\s+class=[\047"]brush:r[\047"]\s*>/../<\/pre>/)' $(TMPFILE) \
		| grep -Ev "<.?pre"


ifdef MAKECMDGOALS
 ifeq ($(MAKECMDGOALS), extractBashCodeFromHtml)
  ifeq ($(HTMLFILE),)
   $(error HTMLFILE variable is undefined. Use: "make extractBashCodeFromHtml HTMLFILE=/path/to/html/file.html")<
  endif
 endif
endif

extractBashCodeFromHtml:
	@perl -ne 'print unless (/<!--/../-->/)' $(HTMLFILE)  > $(TMPFILE)
	@perl -ne 'print if(/<pre\s+class=.*?brush:bash.*?>/../<\/pre>/)' $(TMPFILE)   |  \
	grep -Ev "<.?pre" | perl -npe 's/^\s*//; s/#.*?$$//;' 

