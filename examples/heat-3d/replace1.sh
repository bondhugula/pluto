head -n 114 temp1.c > temp2.c
printf "if(t5%%2==0){" >> temp2.c
sed -n '115,125p' temp1.c >>temp2.c
printf "}else{\n" >> temp2.c
sed -n '115,125p' temp1.c >>temp2.c
printf "}" >> temp2.c
tail -n +126 temp1.c >>temp2.c

sed -in '/if(/,/else{/ s/A\[t5\]/A\[0\]/g' temp2.c
sed -in '/if(/,/else{/ s/A\[t5+1\]/A\[1\]/g' temp2.c
sed -in '/else{/,/icc/ s/A\[t5+1\]/A\[0\]/g' temp2.c
sed -in '/else{/,/icc/ s/A\[t5\]/A\[1\]/g' temp2.c
