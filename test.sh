#!/bin/bash

python -m emigrate.cli construct -i examples/constructor.json -o examples/initial_condition.json
python -m emigrate.cli solve examples/initial_condition.json --output examples/cli_test.hdf5 -t 100.0 -d 1.0
python -m emigrate.cli open examples/cli_test.hdf5 echo -f 10
python -m emigrate.cli open examples/cli_test.hdf5 plot -f 10
