#!/bin/bash
#
# #################
# ## MITObim.sh ###
# #################
#
# Christoph Hahn, May 2012

startiteration=$1
numberofiterations=$2
mode="$3"
strain="$4"
	
MIRArefname="$5"	
MIRAmaf="$6"
MIRAreadpool="$7"        

iteration=$startiteration
let enditeration=$numberofiterations+1     
let seed=$iteration-1

#checking number of input options
if [ $# != 7 ]; then
	echo "Usage: sh test-script.sh [PARAMETERS] [info on initial MIRA mapping assembly]"        
	echo ""	
	echo "Parameters:"
	echo ""	
	echo "[start with iteration] 	<int>		Desired iteration to start with"
	echo "[number of iterations] 	<int>		Total number of iterations desired"
	echo "[mode]			<string>	Mode the wrapper should be used in, SE - Single End, PE - Paired End"
	echo "[strain]			<string>	ID of the strain/population/species under question, e.g. thy1"	
	echo ""
	echo "Required information from MIRA assembly to be used as starting point for iterations"
	echo ""
	echo "[MIRArefname]		<string>	Name of references strain used in Miraproject to be used as starting point for iterations, e.g. der-mt"	
	echo "[PATH to MAF]		<PATH to file>	Full path (absolute) to the Mira *.maf file to be used as starting point for iterations"	
	echo "[PATH to readpool]	<PATH to file>	Full path (absolute) to *.fastq file to be used as readpool in subsequent iterations"
        echo ""	
	echo "NOTE: So far all 7(!) informations are required - the sequence is crucial!!!!!!!!!!!!!!!"
	echo ""
	exit
fi

#checking number of iterations

if [ $startiteration -gt $numberofiterations ]; then
	echo ""
	echo "startiteration has to be smaller than or equal to the total number of iterations - bye!"
	exit

elif [ $startiteration == 1 ]; then
	echo ""
	echo "ok - starting from iteration 1"

elif (($startiteration >= 2)); then
	echo ""
	echo "ok - starting from iteration $startiteration"
else
	echo ""
	echo "somethings wrong with your iterations! - bye!"
	exit
fi


#checking if Mode is set correctly
if [ $mode == "SE" ]; then
	echo ""
	echo "MODE SE - ok"
	uti=no
elif [ $mode == "PE" ]; then
	echo ""
	echo "MODE PE - ok"
	uti=yes
else 
	echo ""
	echo "Select the mode - Singel End (SE) or Paired End (PE"	
	exit
fi

#checking for *.maf file

if [ "$(file $MIRAmaf| awk -F" " '{print $2}')" == "ASCII" ] && [ -f $MIRAmaf ]; then
	echo ""
	echo "*.maf file seems to be in the right format - good!"
else
	echo ""
	echo "something is wrong with the *.maf file - BYE!"
	exit
fi

#checking for readpool file
if [ -f $MIRAreadpool ]; then
	echo ""
        echo "readpool seems not to be a directory - good!"

else
	echo ""
    	echo "readpool file seems to be a directory or something else strange! - BYE!"
	exit
fi



################################## get trimmed readpool ############################################################
	

if [ "$(file $MIRAreadpool| awk -F" " '{print $2}')" == "gzip" ]; then
	echo "copying and unzipping of readpool"
        cp $MIRAreadpool readpool.gz
        gunzip readpool.gz
	MIRAreadpool="../../readpool"
else
        echo "readpool is not zipped - will be used as it is"
         
fi


################################## Starting MITObim ################################################################
#######################################################################################################################

while [  $iteration -lt 2 ]; do
	
		

	echo "  #####################################################   "
	echo "  ########## The iteration is $iteration ##############   "
	echo "  #####################################################   "
	echo "  "

	mkdir Mirabait-iteration$iteration	
	cd Mirabait-iteration$iteration
	mkdir $mode
	cd $mode

########### get new backbone from previous iteration    ########################################

	echo "  "
	echo "  #####################################################   "
	echo "  ########## get new backbone ########################   "
	echo "  #####################################################   "
	echo "  "

	convert_project -f maf -t fasta -A "SOLEXA_SETTINGS -CO:fnicpst=yes" $MIRAmaf $strain-iteration$iteration-seed
	cat $strain-iteration$iteration-seed_$strain.unpadded.fasta | sed 's/@/N/g' > seed.fasta
	
	if [ $(ls -l seed.fasta | awk '{ print $5 }') -gt 0 ] && [ -f seed.fasta ]; then
	
		echo ""
        	echo "seed file has been successfully created!"
		echo ""
	else
		echo ""
    		echo "seed file has not been created or is empty - Stop script!"
		echo "see error message of the convert_project script above for details!"
        	exit
	fi

############################# bait out reads that map onto the backbone #######################################

	echo "  "
	echo "  #####################################################   "
	echo "  ########## MIRABAIT ########################   "
	echo "  #####################################################   "
	echo "  "
		
	mirabait -k 31 -n 1 seed.fasta $MIRAreadpool $strain-iteration$iteration-$mode\_in.solexa

#####checking wheter mirabait yielded output or not
	if [ "$(ls -l  $strain-iteration$iteration-$mode\_in.solexa.fastq | awk '{ print $5 }')" -gt 0 ] && [ -f $strain-iteration$iteration-$mode\_in.solexa.fastq ]; then

		echo ""
       		echo "You caugth some reads - CONGRATULATIONS!"
		echo ""
	else
		echo ""
		echo "baiting reads failed - Stop script!"
		echo "the Mirabait output above might give you an idea what went wrong!"
		echo ""
        	exit
	fi

############################# find mates to mapped reads ######################################################

	if [  $mode == "PE" ]; then		


		echo "  "
    		echo "  #####################################################   "
   		echo "  ########## find mates ########################   "
    		echo "  #####################################################   "
    		echo "  "
		
		mv $strain-iteration$iteration-$mode\_in.solexa.fastq $strain-iteration$iteration-$mode\_in.solexa.fastq-SE
		grep "/1$" $strain-iteration$iteration-$mode\_in.solexa.fastq-SE | sed 's/@HWI/HWI/g' > list1.txt
	   	grep "/2$" $strain-iteration$iteration-$mode\_in.solexa.fastq-SE | sed 's/@HWI/HWI/g' > list2.txt

		sed 's/\/1$/\/2/g' list1.txt > list1-mates.txt
		sed 's/\/2$/\/1/g' list2.txt > list2-mates.txt
		cat list1.txt list2.txt list1-mates.txt list2-mates.txt |sort -n|uniq > $strain-iteration$iteration-$mode.list

		rm list*

		convert_project -f fastq -t fastq -n $strain-iteration$iteration-$mode.list $MIRAreadpool $strain-iteration$iteration-$mode\_in.solexa
	else 
		echo "Mode SE selected - finding pairs omitted!"		
	fi

	rm $strain-iteration$iteration-$mode.list
	if [ "$(ls -l  $strain-iteration$iteration-$mode\_in.solexa.fastq | awk '{ print $5 }')" -gt 0 ] && [ -f $strain-iteration$iteration-$mode\_in.solexa.fastq ]; then

	        echo ""
	        echo "Finding good pairs successfull!"
	        echo ""
        else
              	echo ""
                echo "Error in finding good pairs - Exit!"
                echo ""
                echo ""
                exit
        fi


############################ check number of paired reads #########################################

   	echo "  #####################################################   "
   	echo " 	############ Number of reads in readpool ###########	"
  	echo "  #####################################################   "

   	grep '@HWI' $strain-iteration$iteration-$mode\_in.solexa.fastq |wc -l
	
	echo "  #####################################################   "
  	echo "  "
    	echo "  #####################################################   "


############################ clean ################################################################

	rm $strain-iteration$iteration-seed_AllStrains.*
	rm $strain-iteration$iteration-seed.fasta*
	rm $strain-iteration$iteration-seed_default.*	
	rm $strain-iteration$iteration-seed_$MIRArefname.*
	rm $strain-iteration$iteration-seed_$strain-iteration$seed.padded*
	rm $strain-iteration$iteration-seed_$strain-iteration$seed.unpadded.fasta	


############################ assemble baited reads using backbone build from the assembly of the previous iteration ##################

   	echo "  "
   	echo "  #####################################################   "
 	echo "  ########## MIRA ########################   "
  	echo "  #####################################################   "
    	echo "  "



     	mv seed.fasta $strain-iteration$iteration-$mode\_backbone_in.fasta
	mv $strain-iteration$iteration-seed_$strain-iteration$seed.unpadded.fasta.qual $strain-iteration$iteration-$mode\_backbone_in.fasta.qual

   	mira --project=$strain-iteration$iteration-$mode --job=mapping,genome,accurate,solexa "--noclipping -CL:pec=no" -GE:not=2:mps=4 -MI:somrnl=0 -AS:nop=1 -SB:bsn=seed:bft=fasta:bbq=30 SOLEXA_SETTINGS -CO:msr=no -GE:uti=$uti:tismin=100:tismax=650 -SB:dsn=$strain-iteration$iteration-$mode

################check wheter MIRA ran succesfully	
	if [ $(ls -l $strain-iteration$iteration-$mode\_assembly/$strain-iteration$iteration-$mode\_d_results/$strain-iteration$iteration-$mode\_out.maf | awk '{ print $5 }') -gt 0 ] && [ -f $strain-iteration$iteration-$mode\_assembly/$strain-iteration$iteration-$mode\_d_results/$strain-iteration$iteration-$mode\_out.maf ]; then

		echo ""
        	echo "MIRA ran successfully!"
		echo ""
		echo "Iteration$iteration finished succesfully - CONGRATULATIONS!"
		echo ""
	else
    		echo "MIRA did not run successfully - Stop script! - please refer to the MIRA output above for details on what went wrong!"
        	exit
	fi
	
    	cd ../../

	let iteration=iteration+1
	let seed=seed+1

done


################################## start iterations #########################################################

while [  $iteration -lt $enditeration ]; do
                
	echo "  #####################################################   "
	echo "  ########## The iteration is $iteration ##############   "
	echo "  #####################################################   "
	echo "  "
	echo "  #####################################################   "
	echo "  ########## The seed is $seed ########################   "
	echo "  #####################################################   "
	echo "  "

	mkdir Mirabait-iteration$iteration
                
        cd Mirabait-iteration$iteration
        mkdir $mode
        cd $mode

########### get new backbone from previous iteration    ########################################

	echo "  "
	echo "  #####################################################   "
	echo "  ########## get new backbone ########################   "
	echo "  #####################################################   "
	echo "  "
	
	convert_project -f maf -t fasta -A "SOLEXA_SETTINGS -CO:fnicpst=yes" ../../Mirabait-iteration$seed/$mode\/$strain-iteration$seed-$mode\_assembly/$strain-iteration$seed-$mode\_d_results/$strain-iteration$seed-$mode\_out.maf $strain-iteration$iteration-seed


        cat $strain-iteration$iteration-seed_$strain-iteration$seed-$mode.unpadded.fasta | sed 's/@/N/g' > seed.fasta


	if [ $(ls -l seed.fasta | awk '{ print $5 }') -gt 0 ] && [ -f seed.fasta ]; then

		echo ""
        	echo "seed file has been successfully created!"
		echo ""
	else
		echo ""
    		echo "seed file has not been created or is empty - Stop script!"
		echo "see error message of the convert_project script above for details!"
        	exit
	fi

############################# bait out reads that map onto the backbone #######################################

	echo "  "
	echo "  #####################################################   "
	echo "  ########## MIRABAIT ########################   "
	echo "  #####################################################   "
	echo "  "
		
	mirabait -k 31 -n 1 seed.fasta $MIRAreadpool $strain-iteration$iteration-$mode\_in.solexa

#####checking wheter mirabait yielded output or not
	if [ "$(ls -l  $strain-iteration$iteration-$mode\_in.solexa.fastq | awk '{ print $5 }')" -gt 0 ] && [ -f $strain-iteration$iteration-$mode\_in.solexa.fastq ]; then

		echo ""
		echo "You caugth some reads - CONGRATULATIONS!"
		echo ""
	else
		echo ""
		echo "baiting reads failed - Stop script!"
		echo "the Mirabait output above might give you an idea what went wrong!"
		echo ""
       		exit
	fi

############################# find mates to mapped reads ######################################################

	if [  $mode == "PE" ]; then


                echo "  "
                echo "  #####################################################   "
                echo "  ########## find mates ########################   "
                echo "  #####################################################   "
                echo "  "

                mv $strain-iteration$iteration-$mode\_in.solexa.fastq $strain-iteration$iteration-$mode\_in.solexa.fastq-SE

		grep "/1$" $strain-iteration$iteration-$mode\_in.solexa.fastq-SE | sed 's/@HWI/HWI/g' > list1.txt
                grep "/2$" $strain-iteration$iteration-$mode\_in.solexa.fastq-SE | sed 's/@HWI/HWI/g' > list2.txt

                sed 's/\/1$/\/2/g' list1.txt > list1-mates.txt
                sed 's/\/2$/\/1/g' list2.txt > list2-mates.txt
                cat list1.txt list2.txt list1-mates.txt list2-mates.txt |sort -n|uniq > $strain-iteration$iteration-$mode.list                        

                rm list*

                convert_project -f fastq -t fastq -n $strain-iteration$iteration-$mode.list $MIRAreadpool $strain-iteration$iteration-$mode\_in.solexa
	else
              	echo "Mode SE selected - finding pairs omitted!"

        fi

	rm $strain-iteration$iteration-$mode.list
	
	if [ "$(ls -l  $strain-iteration$iteration-$mode\_in.solexa.fastq | awk '{ print $5 }')" -gt 0 ] && [ -f $strain-iteration$iteration-$mode\_in.solexa.fastq ]; then

                echo ""
                echo "Finding good pairs successfully!"
                echo ""
        else
               	echo ""
                echo "Error in finding good pairs - Exit!"
                echo ""
                echo ""
                exit
        fi


############################ check number of paired reads #########################################

	echo "  #####################################################   "
   	echo " 	############ Number of reads in readpool ###########	"
  	echo "  #####################################################   "

   	grep '@HWI' $strain-iteration$iteration-$mode\_in.solexa.fastq |wc -l
	
	echo "  #####################################################   "
  	echo "  "
    	echo "  #####################################################   "


############################ clean ################################################################

	rm $strain-iteration$iteration-seed_AllStrains.*
	rm $strain-iteration$iteration-seed.fasta*
	rm $strain-iteration$iteration-seed_default.*	
	rm $strain-iteration$iteration-seed_seed.*
	rm $strain-iteration$iteration-seed_$strain-iteration$seed-$mode.padded*
	rm $strain-iteration$iteration-seed_$strain-iteration$seed-$mode.unpadded.fasta


##  	rm $strain-iteration$iteration-$mode_in.solexa.fastq
##  	rm $strain-iteration$iteration-PE.list


############################ assemble baited reads using backbone build from the assembly of the previous iteration ##################

  	echo "  "
  	echo "  #####################################################   "
	echo "  ########## MIRA ########################   "
 	echo "  #####################################################   "
   	echo "  "

    	mv seed.fasta $strain-iteration$iteration-$mode\_backbone_in.fasta
	mv $strain-iteration$iteration-seed_$strain-iteration$seed-$mode.unpadded.fasta.qual $strain-iteration$iteration-$mode\_backbone_in.fasta.qual

	mira --project=$strain-iteration$iteration-$mode --job=mapping,genome,accurate,solexa "--noclipping -CL:pec=no" -GE:not=2:mps=4 -MI:somrnl=0 -AS:nop=1 -SB:bsn=seed:bft=fasta:bbq=30 SOLEXA_SETTINGS -CO:msr=no -GE:uti=$uti:tismin=100:tismax=650 -SB:dsn=$strain-iteration$iteration-$mode

################check wheter MIRA ran succesfully	
	if [ $(ls -l $strain-iteration$iteration-$mode\_assembly/$strain-iteration$iteration-$mode\_d_results/$strain-iteration$iteration-$mode\_out.maf | awk '{ print $5 }') -gt 0 ] && [ -f $strain-iteration$iteration-$mode\_assembly/$strain-iteration$iteration-$mode\_d_results/$strain-iteration$iteration-$mode\_out.maf ]; then

		echo ""
        	echo "MIRA ran successfully!"
		echo ""
		echo "Iteration $iteration finished succesfully - CONGRATULATIONS!"
		echo ""
		else
    		echo "MIRA did not run successfully - Stop script! - please refer to the MIRA output above for details on what went wrong!"
        	exit
	fi

    	cd ../../

	let iteration=iteration+1
	let seed=seed+1
	
done

let finaliteration=$iteration-1

rm readpool

echo ""
echo "YOU HAVE SUCCESSFULLY COMPLETED ITERATION NUMBER $finaliteration, YOU ASKED FOR $numberofiterations, SO THAT'S IT!!!"
echo ""
