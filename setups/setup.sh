#!/bin/bash

# create a copy of the original bashrc file
cp ~/.bashrc ~/.bashrc_backup

# copy the new rc files into your home directory
cp .bashrc ~/.bashrc
cp .inputrc ~/.inputrc
cp .vimrc ~/.vimrc

# activate the settings
bind -f ~/.inputrc
source ~/.bashrc

# copy the data folder into your home directory
mkdir -p ~/data
cp -R -T data ~/data # copy things in data into ~/data (no target directory)

# if [ -d ~/data ]; then
#     echo "~/data/ already exists"
#     cp -R -T data ~/data
# else
#     echo "There's no ~/data/ directory"
#     cp -r data ~
# fi

