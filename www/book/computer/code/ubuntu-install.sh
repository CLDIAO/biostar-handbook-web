# Update your linux package.
# You may need to run this periodically.
sudo apt-get update && sudo apt-get upgrade

# Install a number of required libraries.
sudo apt-get install -y curl build-essential ncurses-dev byacc zlib1g-dev  git cmake

# Install java and python libraries.
sudo apt-get install -y python-pip python-dev libwww-perl libhtml-parser-perl

