VAR=`tr '\n' '-' <tile.sizes  | sed 's/-$//'`
mv temp2.c lb-$VAR.c
mv temp2 lb-ex-$VAR
