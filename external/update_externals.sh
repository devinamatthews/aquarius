#!/bin/bash
for ext in `ls`; do
    if [ -d $ext/.git.bak ]; then
        echo "Updating $ext"
        cd $ext
        mv .git.bak .git
        git pull
        mv .git .git.bak
        cd ..
    fi
done
