#!/bin/bash
set -e

git fetch origin
git merge origin/tmp
git reset HEAD~1
