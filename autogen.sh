#!/bin/bash -e

  BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  export BASE
  
  # ------------------------------------------------
  # (Re)Generate autotools files
  # ------------------------------------------------
  
  echo -e "\n*** Running autotools on pluto ***"
  autoreconf -vi
  
  echo -e "\n*** Running autotools on isl ***"
  (cd isl; ./autogen.sh)
  
  echo -e "\n*** Running autotools on cloog-isl ***"
  (cd cloog-isl; ./autogen.sh)
  
  echo -e "\n*** Running autotools on piplib ***"
  (cd piplib; ./autogen.sh)
  
  echo -e "\n*** Running autotools on polylib ***"
  (cd polylib; ./autogen.sh)
  
  echo -e "\n*** Running autotools on openscop ***"
  (cd openscop; ./autogen.sh)
  
  echo -e "\n*** Running autotools on candl ***"
  (cd candl; ./autogen.sh)
  
  echo -e "\n*** Running autotools on clan ***"
  (cd clan; ./autogen.sh)

