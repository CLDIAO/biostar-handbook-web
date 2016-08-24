
OUTPUT=~/tmp/pyblue/biostar-handbook
WEB=www
REMOTE=www@biostarhandbook.com:/home/www/sites/handbook/www

serve:
	# Default action is to render the pages.
	pyblue -r $(WEB) -p 4000

push:
	git commit -am 'book update'
	git push

html:
	mkdir -p $(OUTPUT)
	pyblue -r $(WEB) -o $(OUTPUT) 

sync: html
	rsync --rsh=ssh -avz $(OUTPUT)/ $(REMOTE)
