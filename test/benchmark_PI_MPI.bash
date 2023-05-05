printf "nb-iter\tnb-process\ttime\n" ;

for nb_iter in $(seq 5 9)
do    
    N=$(python -c "print(10**$nb_iter)");
    # printf "N = $N\n" ;
    for w in $(seq 1 4)
    do
	    # printf "\tNombre processus = $w\n" ;
	    # /usr/bin/time -f "\t\t%E" ./test_pi -np $w --nb-iter $N --quiet;
        printf "$N\t$w\t" ;
        echo $(/usr/bin/time -f "%E" ./test_pi -np $w --nb-iter $N --quiet) -n;
    done
done