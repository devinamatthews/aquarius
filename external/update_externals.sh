#!/bin/bash
for ext in `ls`; do
    if [ -f $ext/.git.bak ]; then
        echo "Updating $ext"
        cd $ext
        tar zxf .git.bak
        git pull
        tar zcf .git.bak .git
        rm -rf .git
        cd ..
    fi
done
