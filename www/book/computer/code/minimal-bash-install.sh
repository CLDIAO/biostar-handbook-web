# Shortcut to the location of the minimal bash.
URL=https://www.biostarhandbook.com/book/computer/minimal-bash.sh

# Always verify the content scripts you get from the internet.
curl $URL | more

# On Mac OSX
curl $URL >> ~/.profile

# On Linux
curl $URL >> ~/.bashrc
