# New version
## 1) Update version in setup.py and commit and push
## 2) Pull request of dev into main
## 3) Make release on GitHub
## 4) Run this code:
rm -r maginator.egg-info/ dist/ build/
python setup.py sdist
twine upload dist/*
