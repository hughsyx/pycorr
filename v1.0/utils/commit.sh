cd ../..
git add `find . \( -name "*.py" -o -name "*.sh" -o -name "*.m" -o -name "*.md" -o -name "Licence*.txt" \)`
git commit -vm ${1}
git push -u origin master
echo "yeah" 