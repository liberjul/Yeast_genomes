grep "T =" ../data/phylo/info/Li_etal2020_strains.txt > ../data/phylo/info/Li_etal2020_type_strains.txt
cat ../data/phylo/info/Li_etal2020_type_strains.txt > ../data/phylo/info/type_strains.txt
cat ../data/phylo/info/extra_type_strains.txt >> ../data/phylo/info/type_strains.txt
python find_types_class.py
