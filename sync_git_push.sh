#!/bin/bash
set -e

# Check that no tmp branch exists in origin
TMP=tmp
if [[ $(git ls-remote --heads origin ${TMP} | wc -l) != "0" ]]
then
    echo "Error: origin/${TMP} already exists"
    exit 1
fi

# Push changes to origin/tmp and delete local changes
git checkout -b tmp
git add .
git commit -m "tmp"
git push origin tmp:tmp
git checkout Dev
git branch -D tmp
