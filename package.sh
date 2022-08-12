rm -r maginator.egg-info/ dist/ build/
python setup.py sdist
python setup.py install
twine upload dist/*
