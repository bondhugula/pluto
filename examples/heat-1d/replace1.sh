head -n 98 temp1.c > temp2.c
printf "if(t3%%2==0){\n" >> temp2.c
sed -n '99,103p' temp1.c >>temp2.c
printf "}else{\n" >> temp2.c
sed -n '99,103p' temp1.c >>temp2.c
printf "}" >> temp2.c
tail -n +104 temp1.c >>temp2.c

sed -in '/if(/,/else{/ s/A\[t3\]/A\[0\]/g' temp2.c
sed -in '/if(/,/else{/ s/A\[t3+1\]/A\[1\]/g' temp2.c
sed -in '/else{/,/icc/ s/A\[t3+1\]/A\[0\]/g' temp2.c
sed -in '/else{/,/icc/ s/A\[t3\]/A\[1\]/g' temp2.c
