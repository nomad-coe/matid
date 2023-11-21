cd ..
python -m pip install --upgrade pip
if [ -f requirements.txt ]; then pip install -r requirements.txt; fi
if [ -f devrequirements.txt ]; then pip install -r devrequirements.txt; fi
