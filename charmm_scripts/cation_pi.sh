#! /bin/csh

set nresidues=`wc cation_pi.candidates|awk '{print $1}'`

touch cation_pi.dat ; rm cation_pi.dat
touch cation_pi.dat2 ; rm cation_pi.dat2

# Generation of a hbonds.dat file
set k=1
while ( $k <= $nresidues ) 

        set resn=`awk -v var=$k 'NR==var {print $1}' cation_pi.candidates|cut -c 1-3`
	set bisresn=`echo $resn|sed 'y/TYRPHE/tyrphe/'`
	if ( $resn != "PHE" && $resn != "TRP" && $resn != "TYR" && $resn != "phe" && $resn != "trp" && $resn != "tyr" ) then
		echo "Wrong residue : must be PHE/TYR/TRP"
		exit
	endif
        set resi=`awk -v var=$k 'NR==var {print $1}' cation_pi.candidates|cut -c 4-6`
 #       set chain=`head -20 ../scr/piplc_protmemb.dyna_desolv.psf |grep CA|awk '{print $8}'|head -1`
        set chain=`head -20 phosphod_bound_to_pure_popc.psf |grep CA|awk '{print $2}'|head -1`
        set skip=1
#        set dd=prod10_wrapped_protmemb.dcd
#        set dd=prod9_wrapped_protmemb.dcd
#        set dd=prod8_wrapped_protmemb.dcd
#        set dd=prod7_wrapped_protmemb.dcd
#        set dd=prod6_wrapped_protmemb.dcd
#        set dd=prod5_wrapped_protmemb.dcd
#        set dd=prod4_wrapped_protmemb.dcd
#        set dd=prod3_wrapped_protmemb.dcd
        set dd=last_50ns.dcd #phospho_d_pure_popc_read.dcd  #
        set dcd=$dd
        set conf_n=`../catdcd4.0/catdcd  -num output/$dcd|grep "Total frames"|awk '{print $NF}'`
	#set inp=`ls ../scr|grep dyna_prod1.cor|awk -F "." '{print $1".dyna_desolv"}'`


	#execution of the charmm input if necessary
	touch ../out/around_cat_pi.$resn$resi.out
#	touch out/around_cat_pi.$resn$resi.out
	set charmmoutputtest=`awk '$1=="NORMAL" && $2=="TERMINATION" {print $0}' ../out/around_cat_pi.$resn$resi.out|wc|awk '{print $1}'`
	if ( $charmmoutputtest != 1 ) then
		echo "Looking what lipids are around residue $resn$resi" 

	        /net/orinoco/apps/charmm/c33b1/exec/gnu/charmm_xxlarge resn:$resn resi:$resi skip:$skip conf_n:$conf_n <around_cat_pi.inp> ../out/around_cat_pi.$resn$resi.out
	else
		echo "...Output out/around_cat_pi.$resn$resi.out already finished"
	endif




        awk '$8=="MEMB" {print $0}' ../out/around_cat_pi.$resn$resi.out|sort -k 2 > $resn$resi.close.tmp

        touch $resn$resi.close.dat ; rm $resn$resi.close.dat
        while ( `wc $resn$resi.close.tmp|awk '{print $1}'` != 0 )
                head -1 $resn$resi.close.tmp >>  $resn$resi.close.dat
                set a=`tail -1 $resn$resi.close.dat|awk '{print $1}'`
                awk -v vara=$a '$1!=vara {print $0}' $resn$resi.close.tmp>$resn$resi.close.tmp2
                mv $resn$resi.close.tmp2 $resn$resi.close.tmp
        end
        rm $resn$resi.close.tmp




        set ncases=`wc $resn$resi.close.dat|awk '{print $1}'`
        set j=1
        while ( $j <= $ncases )
                set lipid=`awk -v var=$j 'NR==var {print $9}'  $resn$resi.close.dat`
                set lipidname=`awk -v var=$j 'NR==var {print $3}'  $resn$resi.close.dat`

		touch $resn$resi.$lipidname$lipid.yesno ; rm $resn$resi.$lipidname$lipid.yesno

		if ( $1 != "nots" ) then	
			echo "Computing timeseries between carbons of $resn$resi and N+ of $lipidname$lipid"
	#		sed 's/XXXDCD/'$dcd'/' ../lib/timeseries.tmp > timeseries.$resn$resi.$lipid.inp
			sed 's/XXXDCD/'$dcd'/' timeseries.tmp > timeseries.$resn$resi.$lipid.inp

			touch $bisresn$resi.$lipid.cationpi.ts
			set charmmoutputtest=`wc $bisresn$resi.$lipid.cationpi.ts|awk '{print $1}'`
			if ( $charmmoutputtest < $conf_n ) then
				#/net/orinoco/apps/charmm/c33b1/exec/gnu/charmm_xxlarge  type:cationpi skip:$skip resn:$resn resi:$resi chain:PROA lipid:$lipid stop:$conf_n < timeseries.$resn$resi.$lipid.inp > ../out/timeseries.$resn$resi.$lipid.out
				/net/orinoco/apps/charmm/c38b2/exec/gnu_xxlarge/charmm_xxlarge type:cationpi skip:1 resn:$resn resi:$resi chain:PROA lipid:$lipid stop:$conf_n < timeseries.$resn$resi.$lipid.inp > ../out/timeseries.$resn$resi.$lipid.out
			else
				echo "...Job already finished"
			endif
		endif

		# first criteria : all carbons < 7A
		awk '{if ($2<=7 && $3<=7 && $4<=7 && $5<=7 && $6<=7 && $7<=7 && $8<=7 && $9<=7 && $10<=7 && NR>14) {print $0} else {print 0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10}}' $bisresn$resi.$lipid.cationpi.ts > $resn$resi.$lipidname$lipid.tmp
		#####awk -v var=$conf_n '$2<=7 && $3<=7 && $4<=7 && $5<=7 && $6<=7 && $7<=7 && $8<=7 && $9<=7 && $10<=7 && NR>14 {print $0}' *$resi.$lipid.cationpi.ts > $resn$resi.$lipidname$lipid.tmp

		set nlines=`wc $resn$resi.$lipidname$lipid.tmp|awk '{print $1}'`
		# second criteria : max-min<1.5
		if ( $resn == "PHE" || $resn == "TYR" ) then
                        set l=1
                        set m=0
                        while ( $l <= $nlines )
                                set mini=`awk -v var=$l 'NR==var {printf("%5.2f\n%5.2f\n%5.2f\n%5.2f\n%5.2f\n%5.2f\n",$2,$3,$4,$5,$6,$7)}' $resn$resi.$lipidname$lipid.tmp|sort --numeric-sort -k 1|head -1`
                                set diff=`awk -v var=$l 'NR==var {printf("%5.2f\n%5.2f\n%5.2f\n%5.2f\n%5.2f\n%5.2f\n",$2,$3,$4,$5,$6,$7)}' $resn$resi.$lipidname$lipid.tmp|sort --numeric-sort -k 1|tail -1|awk -v var=$mini '{printf("%5d\n",($1-var)*1000)}'`
                                if ( $diff <= 1500 ) then
					echo "1" >> $resn$resi.$lipidname$lipid.yesno
                                        @ m = ( $m + 1 )
				else
                                        echo "0" >> $resn$resi.$lipidname$lipid.yesno
                                endif
                                @ l = ( $l + 1 )
                        end
			goto finish
		endif

		if ( $resn == "TRP" ) then
			set l=1
			set m=0
			while ( $l <= $nlines )
				set mini=`awk -v var=$l 'NR==var {printf("%5.2f\n%5.2f\n%5.2f\n%5.2f\n%5.2f\n%5.2f\n%5.2f\n%5.2f\n%5.2f\n",$2,$3,$4,$5,$6,$7,$8,$9,$10)}' $resn$resi.$lipidname$lipid.tmp|sort --numeric-sort -k 1|head -1`
                                set diff=`awk -v var=$l 'NR==var {printf("%5.2f\n%5.2f\n%5.2f\n%5.2f\n%5.2f\n%5.2f\n%5.2f\n%5.2f\n%5.2f\n",$2,$3,$4,$5,$6,$7,$8,$9,$10)}' $resn$resi.$lipidname$lipid.tmp|sort --numeric-sort -k 1|tail -1|awk -v var=$mini '{printf("%5d\n",($1-var)*1000)}'`
				if ( $diff <= 1500 ) then
                                        echo "1" >> $resn$resi.$lipidname$lipid.yesno
					@ m = ( $m + 1 )
                                else
                                        echo "0" >> $resn$resi.$lipidname$lipid.yesno
				endif
				@ l = ( $l + 1 )
			end
		endif

                echo $resn$resi $lipid $lipidname $nlines $m

		finish:
		set cationpicountb=`echo $m|awk -v var=$conf_n '{{printf("%4d",$1/var*100)}}'`

		echo " > $resn$resi / $lipidname$lipid : $cationpicountb %" >> cation_pi.dat

                @ j = ( $j +  1 )
        end

        paste $resn$resi*yesno|awk '{print NR" "$1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20}' > $resn$resi.yesno.all
        set avea=`awk -v var=$conf_n 'BEGIN {sumdiff=0} {sumdiff+=$2} END {printf("%4.1f\n",sumdiff/var*100)}' $resn$resi.yesno.all`
	awk '{if ($2>0) {print $1" "$2" "1" "$2-1} else {print $1" "$2" "0" "0}}' $resn$resi.yesno.all > proutetmp
	mv proutetmp $resn$resi.yesno.all
        set aveb=`awk -v var=$conf_n 'BEGIN {sumdiff=0} {sumdiff+=$3} END {printf("%4.1f\n",sumdiff/var*100)}' $resn$resi.yesno.all`
        set avec=`awk -v var=$conf_n 'BEGIN {sumdiff=0} {sumdiff+=$4} END {printf("%4.1f\n",sumdiff/var*100)}' $resn$resi.yesno.all`

	echo " > $resn$resi / all : $avea $aveb $avec" >> cation_pi.dat2
	echo ""
        echo " > $resn$resi / all : $avea $aveb $avec"
	echo ""

	echo "-100 100" > $resn$resi.yesno.all.downsized
	set r=100
	while ( $r <= $conf_n )
		awk -v var=$r 'NR>(14+var-100) && NR<=(14+var) {print $3}' $resn$resi.yesno.all>$resn$resi.yesno.all.caca
		set testcaca=`grep 1 $resn$resi.yesno.all.caca|wc|awk '{print $1}'`
		if ( $testcaca >= 50 ) then
			echo "$r 1" >> $resn$resi.yesno.all.downsized
		endif
		#set bgbg=`~/.scripts/average_stdev.sh 3 $resn$resi.yesno.all.caca|awk '{printf("%2d\n",$1)}'`
		#echo "$r $bgbg" >> $resn$resi.yesno.all.downsized
		@ r = ( $r + 100 )
	end
	rm $resn$resi.yesno.all.caca

	@ k = ( $k + 1 )
end
