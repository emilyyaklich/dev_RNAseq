#!/bin/bash
#SBATCH --job-name=write_interpro_submission
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=14
#SBATCH --mem=5gb
#SBATCH --time=48:00:00
#SBATCH --output=write_interpro_%j.out
#SBATCH --error=write_interpro_%j.error
#SBATCH --mail-user=ely67071@uga.edu
#SBATCH --mail-type=ALL


work_dir="/scratch/ely67071/sunflower_inflo_dev_data_b3/interpro"
cd ${work_dir}

for run_number in *.fa

do
    

        echo "#!/bin/bash" > interpro_scripts/${run_number}.sh

        echo "#SBATCH --partition=batch" >> interpro_scripts/${run_number}.sh

        echo "#SBATCH --nodes=1" >> interpro_scripts/${run_number}.sh

        echo "#SBATCH --ntasks=16" >> interpro_scripts/${run_number}.sh

        echo "#SBATCH --mem=50G" >> interpro_scripts/${run_number}.sh

        echo "#SBATCH --time=48:00:00" >> interpro_scripts/${run_number}.sh

        echo "#SBATCH --output=%x_%j.out" >> interpro_scripts/${run_number}.sh

        echo "#SBATCH --error=%x_%j.err" >> interpro_scripts/${run_number}.sh

        echo "#SBATCH --mail-type=END,FAIL" >> interpro_scripts/${run_number}.sh

        echo "#SBATCH --mail-user=ely67071@uga.edu" >> interpro_scripts/${run_number}.sh

        echo "">> interpro_scripts/${run_number}.sh

        echo "ml InterProScan/5.64-96.0-foss-2022b" >> interpro_scripts/${run_number}.sh

        echo "sh interproscan.sh -goterms -i ${run_number}" >> interpro_scripts/${run_number}.sh       

        chmod 755 interpro_scripts/${run_number}.sh

done    

 
cp *.fa interpro_scripts/
