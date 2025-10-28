for t in 1 5 10 11 12 13 14 15 16 17 18 19 20 30 40 50
do
    scaled_t=$(awk "BEGIN {print $t * 3.575}")
    gtime -f "%E %M" ./prepro_verification_aids -f FORI -w . -l 1800 -t 12 -r 112358 -s 10 -v ${scaled_t} --flat 1 -p 1 
done
for t in 1 5 10 11 12 13 14 15 16 17 18 19 20 30 40 50
do
scaled_t=$(awk "BEGIN {print $t * 3.575}")
gtime -f "%E %M" ./prepro_verification_muta -f FORI -w . -l 1800 -t 12 -r 112358 -s 10 -v ${scaled_t} --flat 1 -p 1
done
for t in 1 5 10 11 12 13 14 15 16 17 18 19 20 30 40 50
do
scaled_t=$(awk "BEGIN {print $t * 8.375}")
gtime -f "%E %M" ./prepro_verification_prot -f FORI -w . -l 1800 -t 12 -r 112358 -s 10 -v ${scaled_t} --flat 1 
done
