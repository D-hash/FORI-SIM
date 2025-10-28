for t in 1 5 10 11 12 13 14 15 16 17 18 19 20 30 40 50
do
gtime -f "%E %M" ./prepro_verification_aids -f FORI -w . -l 1800 -t 12 -r 112358 -s 10 -v ${t} --flat 1 -p 1 -u 1
done
for t in 1 5 10 11 12 13 14 15 16 17 18 19 20 30 40 50
do
gtime -f "%E %M" ./prepro_verification_muta -f FORI -w . -l 1800 -t 12 -r 112358 -s 10 -v ${t} --flat 1 -p 1 -u 1
done
