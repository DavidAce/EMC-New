#!/bin/bash

# trap ctrl-c and call ctrl_c()
trap ctrl_c INT

function ctrl_c() {
        echo "  **  Trapped CTRL-C"
}

./EMC xy_old.dat xy_new.dat boundaries.dat
echo "Finished Job"