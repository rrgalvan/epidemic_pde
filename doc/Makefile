TESTS=tests_description
PDFVIEW=atril

$(TESTS).pdf: $(TESTS).md
	pandoc $(TESTS).md -o $(TESTS).pdf

view: $(TESTS).pdf
	$(PDFVIEW) $(TESTS).pdf

clean:
	rm -f $(TESTS).pdf *.aux *.out *.fls
