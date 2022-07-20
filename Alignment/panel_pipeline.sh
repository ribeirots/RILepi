# pipeline to obtain panel files for both RIL sets
# tribeiro@wisc.edu


for c in X 2L 2R 3L 3R
do
	for l in ZI418N ZI251N EF43N FR320N
	do
		mv ${l}_Chr${c}_diploid.fasta ${l}_Chr${c}.fas
	done
done

mkdir FR_panel
mkdir EF_panel

for c in X 2L 2R 3L 3R
do
	perl find_snps.pl FR320N_Chr${c}.fas ZI418N_Chr${c}.fas > FR_panel/${c}.panel
	perl find_snps.pl EF43N_Chr${c}.fas ZI251N_Chr${c}.fas > EF_panel/${c}.panel
done
