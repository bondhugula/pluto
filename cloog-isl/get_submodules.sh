#!/bin/sh
git submodule init
git submodule update
(cd isl; git submodule init; git submodule update)
