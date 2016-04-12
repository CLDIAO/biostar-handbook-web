
OUTPUT=~/tmp/pyblue/biostar-handbook
WEB=www
REMOTE=www@biostarhandbook.com:/home/www/sites/handbook

serve:
	# Default action is to render the pages.
	pyblue -v -r $(WEB)

html:
	mkdir -p $(OUTPUT)
	pyblue -r $(WEB) -o $(OUTPUT)

web: html
	rsync --rsh=ssh -avz $(OUTPUT)/ $(REMOTE)