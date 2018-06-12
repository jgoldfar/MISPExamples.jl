SHELL=/bin/bash
PWD:=$(shell pwd)

LATEX=latexmk -bibtex -pdf

TEX_SOURCES=$(wildcard *.tex)
TEX_SOURCES_BASENAME=$(basename $(TEX_SOURCES))
ROOTCOMPILE=main

CHKTEXARGS=-n3
CHKTEX=chktex -q $(CHKTEXARGS)
UNAME=$(shell uname -s)

all: main.pdf

# Visualize output using Skim or Evince depending on platform
ifeq ($(UNAME), Darwin)
viz-%: %.pdf
	/usr/bin/osascript -e "set theFile to POSIX file \"$(PWD)/$<\" as alias" -e "set thePath to POSIX path of theFile" -e "tell application \"Skim\"" -e "  set theDocs to get documents whose path is thePath" -e "  try" -e "    if (count of theDocs) > 0 then revert theDocs" -e "  end try" -e "  open theFile" -e "end tell" &
endif
ifeq ($(UNAME), Linux)
viz-%: %.pdf
	evince "$(PWD)/$<" &
endif

#TODO: Re-add spellcheck-all
check: check-srcs

clean: clean-srcs clean-fmt clean-check clean-spell

clean-all: clean clean-all-srcs
	$(RM) main-flat.tex

.PHONY: list
list:
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' | xargs

## Dependencies
numerics-model-1-deps=$(addprefix verif/TSE-,nu.tex g.tex p.tex)

verif/TSE-%.tex: verif/Makefile
	$(MAKE) -C $(dir $<) $(notdir $@)

.SECONDARY: $(numerics-model-1-deps)

#TODO: Consider building dependency list using latexmk -deps-out=file
## Generate PDF files (and corresponding logfiles) from TeX files
%.pdf %.log: %.tex $(filter-out %.tex,$(TEX_SOURCES)) $(numerics-model-1-deps)
	$(LATEX) $<

## Clean targets for TeX files
clean-src-%: %.tex
	$(LATEX) -c $<

clean-srcs: $(addprefix clean-src-,$(TEX_SOURCES_BASENAME))

clean-all-src-%: %.tex
	$(LATEX) -C $<

clean-all-srcs: $(addprefix clean-all-src-,$(TEX_SOURCES_BASENAME))

## Automatic indentation
%.bak: %.tex
	echo "Indenting $<"
	latexindent -w -l $<

fmt: $(patsubst %.tex,%.bak,$(filter-out main-full.tex,$(TEX_SOURCES)))

clean-fmt:
	$(RM) *.bak
	$(RM) indent.log

## Check targets for various issues we can ascertain from the logfile
check-logs: $(addprefix check-log-,$(ROOTCOMPILE))

check-log-%: %.log
	@grep "Overfull" $*.log | grep "hbox" > $*.overfull.log || echo "No overfull hboxes."
	test ! -s "$*.overfull.log"
	@$(RM) $*.overfull.log
	@grep "LaTeX Warning: There were undefined references." $*.log > $*.undef.log || echo "No undefined references."
	test ! -s "$*.undef.log"
	@$(RM) $*.undef.log

## Check target using linter
lint-%.out: %.tex
	$(CHKTEX) $< 2>/dev/null | tee $@

check-src-%: lint-%.out
	test ! -s $<

check-srcs: $(addprefix check-src-,$(TEX_SOURCES_BASENAME))

clean-check-%: %.tex
	$(RM) lint-$*.out

clean-check: $(addprefix clean-check-,$(TEX_SOURCES_BASENAME))

# TODO: Move all configuration variables like this one to a section above.
## Target to compile other (non-main) files (intended for CI only)
l-check-compile: $(addsuffix .pdf,$(ROOTCOMPILE))

## Spellcheck target
spell/%.log: %.tex
	< $< aspell -W 2 --personal=./.aspell.en.pws \
				--add-tex-command="qquad op" \
				--add-tex-command="eqref p" \
				--add-tex-command="cref p" \
				--add-tex-command="Cref p" \
				--add-tex-command="cite p" \
				--add-tex-command="cpageref p" \
				--add-tex-command="citep p" \
				--add-tex-command="todo op" \
				--add-tex-command="spt o" \
				--add-tex-command="emph P" \
				--add-tex-command="lnorm p" \
				--add-tex-command="qedhere p" \
				--add-tex-command="label p" \
				--add-tex-command="frac pp" \
				--add-tex-command="leq o" \
				--add-tex-command="text p" \
				--add-tex-command="numri o" \
				--dont-tex-check-comments \
				--dont-suggest --lang=en --mode=tex -a &> $@

spellcheck-all: $(addprefix spell/,$(patsubst %.tex,%.log,$(TEX_SOURCES)))

clean-spell:
	$(RM) spell/*.log

check: check-srcs check-logs

## Generate inlined files using custom latexpand script
%-flat.tex: %.tex latexpand $(filter-out %-flat.tex,$(SOURCES))
	./latexpand --strip-comments $< > $@
