#!/bin/bash
. assert.sh
assert_raises 'emigrate construct -i examples/constructor.json -o examples/initial_condition.json'
assert_raises 'emigrate load examples/initial_condition.json solve --output examples/cli_test.hdf5 -t 10.0 -d 1.0'
assert_raises 'emigrate load examples/cli_test.hdf5 echo -f 5'
assert_raises 'emigrate load examples/cli_test.hdf5 plot -f 5 examples/test_plot.png'
emigrate load examples/cli_test.hdf5 echo -f 5
assert_end
