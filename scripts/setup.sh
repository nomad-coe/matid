cd ..
python -m pip install --upgrade pip
if [ -f devrequirements.txt ]; then pip install -r devrequirements.txt; fi
pip install .
