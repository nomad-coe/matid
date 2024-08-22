cd ..
rm -rf dist
pipx run build --sdist
# Before uploading with twine, fetch the wheels from github action arfifact and
# put them inside the 'dist' folder
# pipx run twine upload dist/*
