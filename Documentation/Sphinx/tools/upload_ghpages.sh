#!/bin/bash - 

repo="https://github.com/DiffusionMRITool/DiffusionMRITool.github.io.git"

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
HTML_DIR=${DIR}/../_build/html

if [ ! -e "$HTML_DIR/index.html" ]; then
  echo "$HTML_DIR does not contain an index.html"
  exit 1
fi

cd $HTML_DIR
rm -rf .git
git init
git add * -f
touch .nojekyll
git add .nojekyll
git commit -a -m "dmritool homepage by sphinx. No history."
git remote add origin $repo
git push origin master --force
rm -rf .git 
