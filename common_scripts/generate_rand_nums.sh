for i in {1..100};
do
   prev=(0 50)
   loop=1
   while [ $loop -lt 12 ];
   do
      num=$(($RANDOM%22+1))
	  num_in_prev=0
	  for p in ${prev[@]}
      do
         if [ $p -eq $num ] ; then
            num_in_prev=1
         fi
      done
	  if [ $num_in_prev -eq 0 ] ; then
	     prev+=("$num")
		 echo $num >> /data/eichertd/som_vn_data/chroms_rand_$i
		 loop=$[$loop+1]
	  fi
   done
   oppchroms=$(comm -3 <( echo ${prev[@]} | tr " " "\n" | sort ) <( echo "0 50 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22" | tr " " "\n" | sort ))
   for j in $oppchroms;
   do
      echo $j >> /data/eichertd/som_vn_data/notchroms_rand_$i
   done
done 