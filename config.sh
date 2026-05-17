S1="sample1"
S2="sample2"

MC1="data/sample1.mcool"
MC2="data/sample2.mcool"

RES=10000

URI1="${MC1}::/resolutions/${RES}"
URI2="${MC2}::/resolutions/${RES}"

CHR="chr2"
START=10000000
END=11000000
REGION="${CHR}:${START}-${END}"

THREADS=4

WINDOWS="50000 100000 150000 200000 300000"
MAIN_WINDOW=100000
