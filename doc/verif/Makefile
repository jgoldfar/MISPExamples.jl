MAXIMACMD = maxima --very-quiet

StandardMaximaFix = sed -f fix-maxima-tex-out-1.sed | tr -d '\n'\
	| sed -f fix-maxima-tex-out-nothinspace.sed \
	| sed -f fix-maxima-tex-otherspacing.sed \
	| sed 's/ ,/,/g'

maxima-TSE-%.tex: TikhonovSamarskii.mac
	$(MAXIMACMD) --batch-string 'batch("$<"); load("mactex-utilities"); with_stdout("$@", print(tex1($*)));'$

TSE-%.tex: maxima-TSE-%.tex
	cat $< | $(StandardMaximaFix) > $@

clean-TSE:
	$(RM) TSE-*.tex

clean: clean-TSE
