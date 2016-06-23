#!/bin/bash - 

repo="git@github.com:DiffusionMRITool/dmritool-doxygen.git"

USAGE="$0 DMRITOOL_BUILD_DIR"
BUILD_DIR=$1
if [ -z "$BUILD_DIR" ]; then
 echo $USAGE
 exit 1
fi

HTML_DIR=${BUILD_DIR}/Documentation/Doxygen/html

if [ ! -e "$HTML_DIR/index.html" ]; then
  echo "$HTML_DIR does not contain an index.html"
  exit 1
fi

cd $HTML_DIR
rm -rf .git
git init
git checkout -b gh-pages
git add * -f
touch .nojekyll
git add .nojekyll
git commit -a -m "Doxygen document. No history."
git remote add origin $repo
git push origin gh-pages --force
rm -rf .git 
