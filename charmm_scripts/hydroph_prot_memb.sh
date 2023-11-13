#! /bin/csh

set charmm="/net/orinoco/apps/charmm/c33b1/exec/gnu/charmm_xxlarge"
set chain=`head phosphod_bound_to_pure_popc.crd|grep CA|awk '{print $8}'|head -1|sed 'y/PRO/pro/'`

if ( $1 == "" ) then
	echo "Give the number of the DCDfile in first argument"
	exit
endif
if ( $2 != "" ) then
        goto sum
endif

#set dcd="prod${1}_wrapped.dcd"
#set dcd="51_100_ns_wrapped.dcd"
#set dcd="prod3_wrapped.dcd"
#set dcd="prod4_wrapped.dcd"
#set dcd="prod5_wrapped.dcd"
#set dcd="prod6_wrapped.dcd"
#set dcd="prod7_wrapped.dcd"
#set dcd="prod8_wrapped.dcd"
#set dcd="prod9_wrapped.dcd"
#set dcd="prod10_wrapped.dcd"
set dcd="phospho_d_pure_popc_read.dcd"

if ( -f "output/$dcd" ) then
	set conf_n=`../catdcd4.0/catdcd -num output/$dcd|grep "Total frames:"|awk '{print $NF}'`
	if ( $1 == 1 ) then
                set confbn=`../catdcd4.0/catdcd -num output/$dcd|grep "Total frames:"|awk '{print $NF-1000}'`
	else
		set confbn=`../catdcd4.0/catdcd -num output/$dcd|grep "Total frames:"|awk '{print $NF}'`
	endif
else
	echo "output/$dcd : File not found"
	exit
endif

echo "####################"
echo "# $dcd "
echo "####################"
# PRO7 43   ILE  HD3  -  MEMB 71   DMPC H5S                0.10      0.24

touch BIGAVERAGE_hydr.$dcd;rm BIGAVERAGE_hydr.$dcd
foreach case ( `cat hydroph_candidates` )

        set resn=`echo $case|cut -c 1-3`
        set resi=`echo $case|cut -c 4-6`

        set output="../out/hydroph_prot_memb.$dcd.$resn$resi.out"

	echo "`date` : $resn$resi"

	touch $output
	set charmmoutputtest=`awk '$1=="NORMAL" && $2=="TERMINATION" {print $0}' $output|wc|awk '{print $1}'`
        if ( $charmmoutputtest != 1 ) then
		$charmm resn:$resn resi:$resi conf_n:$conf_n dcd:$dcd < hydroph_prot_memb.inp > $output
	endif

        set time=`sed '1,/RUNNING FROM STEP/d' $output|grep "Resolution"|awk '{print $3}'|sed 's/ps.//'|awk '{print $1*100}'`
        set occupancy=`tail -50 $output|grep Total|grep frame|awk '{printf("%4.2f\n",$3)}'`
        if ( $occupancy == "" ) then
                set occupancy="0.00"
        endif

        foreach output ( `ls ../out/hydroph_prot_memb.${dcd}.$resn$resi.out` )
	        set time=`sed '1,/RUNNING FROM STEP/d' $output|grep "Resolution"|awk '{print $3}'|sed 's/ps.//'|awk '{print $1*100}'`

                awk 'NF==11 {print $0}' $output|grep MEMB|awk '{print $1" "$2" "$3" "$4" "$7" "$8" "$9" "($11-$10)*100" "$11*100}'|sort -n -k 8 > $output.dat1

                set testsize=`ls -l $output.dat1|awk '{print $5}'`
                if ( $testsize == 0 ) then
                        set lifetime="0.00"
                        goto zapthat
                else
                        set nlines=`wc $output.dat1|awk '{print $1+1}'`
                        if ( $nlines == 2 ) then
                                set lifetime=`head -1 $output.dat1|awk -v var=$time '{printf("%5.2f\n",($9-$8+1)/var)}'`
                                goto zapthat
                endif

                set i=1
                set j=2
                touch $output.dat2 ; rm $output.dat2
                set end=0
                set final=`tail -50 $output.dat1|sort -n -k 9|tail -1|awk '{print $NF-1}'`
                while ( $end <= $final )
			# PRO7 43   ILE  HD3  -  MEMB 71   DMPC H5S                0.10      0.24
                        set begin=`awk -v var=$i 'NR==var {print $(NF-1)}' $output.dat1`
                        set end=$begin

                        loop:

                        set b=`awk -v var=$i 'NR==var {print $(NF)}' $output.dat1`
                        set c=`awk -v var=$j 'NR==var {print $(NF-1)}' $output.dat1`
                        @ j = ( $j + 1 )
                        @ i = ( $i + 1 )

                        if ( $c == "" ) then
                                if ( $b > $end ) then
                                        set end=$b
                                        goto truc
                                else
                                        goto truc
                                endif
                        endif

                        if ( $b >= $c ) then
                                if ( $b <= $end ) then
                                        goto loop
                                else
                                        set end=$b
                                        goto loop
                                endif
                        else
                                if ( $c <= $end ) then
                                        goto loop
                                endif
                                if ( $b > $end ) then
                                        set end=$b
                                endif
                        endif
                        truc:
                        set duration=`echo $begin|awk -v var=$end '{print var-$1+1}'`
                        echo "$begin $end $duration $i $j" >> $output.dat2
                end

                set lifetime=`awk -v var=$time 'BEGIN {sum1=0} {sum1+=$3} END {printf("%5.2f\n",sum1/var*100)}' $output.dat2`

                zapthat:

                echo "$resn$resi typea typeb $lifetime $time $occupancy" >> BIGAVERAGE_hydr.$dcd
		echo "         -> $lifetime % / $occupancy"
	end

end

exit

sum:

set divid=`paste BIGAVERAGE*dcd|head -1|awk '{printf("%10d\n",$5+$11+$17+$23+$29+$35+$41+$47+$53+$59+$65+$71+$77+$83+$89+$95+$101+$107+$113+$119+$125+$131+$137+$143+$149)}'`
paste BIGAVERAGE*dcd|grep -v -y "val[1-9]* sc"|grep -v -y "leu[1-9]* sc"|grep -v -y "gly[1-9]* sc"|grep -v -y "pro[1-9]* sc"|grep -v -y "met[1-9]* sc"|grep -v -y "ala[1-9]* sc"|awk -v var=$divid '{print $1" "$2" "$3" "($4*$5+$10*$11+$16*$17+$22*$23+$28*$29+$34*$35+$40*$41+$46*$47+$52*$53+$58*$59+$64*$65+$70*$71+$76*$77+$82*$83+$88*$89+$94*$95+$100*$101+$106*$107+$112*$113+$118*$119+$124*$125+$130*$131+$136*$137+$142*$143+$148*$149+$154*$155+$160*$161+$166*$167)/var" "($6*$5+$12*$11+$18*$17+$24*$23+$30*$29+$36*$35+$42*$41+$48*$47+$54*$53+$60*$59+$66*$65+$72*$71+$78*$77+$84*$83+$90*$89+$96*$95+$102*$101+$108*$107+$114*$113+$120*$119+$126*$125+$132*$131+$138*$137+$144*$143+$150*$149+$156*$155+$162*$161+$168*$167)/var}'|sort -n -k 4|awk '$4!=0 {print $0}' > BIGAVERAGE_all

