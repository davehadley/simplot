#!/bin/bash

TESTPACKAGE=${SIMPLOTROOT}/test
python -m unittest discover ${TESTPACKAGE}
