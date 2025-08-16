# download the all core-nt database from NCBI
for i in $(seq -f "%02g" 36 74); do
	echo "$i"
	wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/core_nt.${i}.tar.gz
	tar -xzvf core_nt.${i}.tar.gz
done

