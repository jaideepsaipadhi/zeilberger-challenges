# Independent check of the 10 x 42 triple

    gcc -O3 -march=native -o ichomp ichomp.c
    ./ichomp 10 42 10 42        # ~2.5-3 h on one core, needs ~2 GB RAM

Expected output:

    10 x 42 winning bites (row,col, 1-based): (5,36) (7,30) (8,26)   total 3

Census of every bar inside the box (non-unique bars only):

    ./ichomp 10 42 > census_10x42.txt
