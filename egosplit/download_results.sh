#!/usr/bin/env bash

cd results_compute
scp -r -C awiebigke@compute5.iti.kit.edu:Code/networkit/egosplit/results/ .
mv results/*.result .
rm -r results