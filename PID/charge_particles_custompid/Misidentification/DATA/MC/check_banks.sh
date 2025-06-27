for dir in /lustre24/expphy/volatile/clas12/dmat/test/febCarbonboth/*; do
    if [ -f "$dir/sidis_mc-master/r_ttest3S.hipo" ]; then
        echo "Checking banks in $dir/sidis_mc-master/r_ttest3S.hipo"
        /u/scigroup/cvmfs/hallb/clas12/sw/noarch/coatjava/11.1.1/bin/hipo-utils -banks "$dir/sidis_mc-master/r_ttest3S.hipo" | grep "MC::RecMatch"
    fi
done
