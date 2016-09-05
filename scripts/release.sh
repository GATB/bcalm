release=$(cat ../VERSION)
git tag -a $release -f
git push --tags -f

