head -n 99 temp1.c > temp2.c
printf "if(t4%%2==0){" >> temp2.c
sed -n '100,108p' temp1.c >>temp2.c
printf "}else{\n" >> temp2.c
sed -n '100,108p' temp1.c >>temp2.c
printf "}" >> temp2.c
tail -n +109 temp1.c >>temp2.c

sed -in '/if(/,/else{/ s/A\[t4\]/A\[0\]/g' temp2.c
sed -in '/if(/,/else{/ s/A\[t4+1\]/A\[1\]/g' temp2.c
sed -in '/else{/,/icc/ s/A\[t4+1\]/A\[0\]/g' temp2.c
sed -in '/else{/,/icc/ s/A\[t4\]/A\[1\]/g' temp2.c
