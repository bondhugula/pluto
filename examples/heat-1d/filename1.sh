VAR=`tr '\n' '-' <tile.sizes  | sed 's/-$//'`
mv temp2.c tp-$VAR.c
mv temp2 par
