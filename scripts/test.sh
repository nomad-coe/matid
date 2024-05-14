#!/usr/bin/env bash
cd ../tests
export COVERAGE_FILE="../reports/coverage/.coverage"
coverage run --source="matid" testrunner.py
unittest=$?
coverage run -m --source="matid" --append pytest
pytest=$?
if [ "$unittest" != 0 ] || [ "$pytest" != 0 ]; then
    exit 1
fi
cd ../reports/coverage
export COVERAGE_FILE=".coverage"
coverage json -o coverage.json
