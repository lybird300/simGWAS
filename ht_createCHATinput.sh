#!/bin/sh

index=0
for cv in {1..20..1} # {start..end..increment}
do
  for cc in 0 1
    do
      for grr in 0 1
        do
          ((index++))

          echo "#!/bin/sh" >> crtCHATinput${index}.sh

          echo "#SBATCH -N 1" >> crtCHATinput${index}.sh
         
          echo "#SBATCH -n 1" >> crtCHATinput${index}.sh
          
          echo "#SBATCH --job-name=simGWA" >> crtCHATinput${index}.sh
          
          echo "#SBATCH --mem-per-cpu=4000" >> crtCHATinput${index}.sh       

          echo "cd /projects/sequence_analysis/vol4/linly/workspace/simGWASNew/" >> crtCHATinput${index}.sh

          echo "java -jar simGWAStag.jar simGWAS.properties ${cv} ${cc} ${grr}" >> crtCHATinput${index}.sh
          
          #echo "sleep 2" >> crtCHATinput${index}.sh
          
          # After editing a shell script, change permission to make it executable by
          # chmod +x script_name
          # Then we can run the script simply by
          # ./script_name
          chmod +x crtCHATinput${index}.sh
        done
    done
done     
exit 0                                  


