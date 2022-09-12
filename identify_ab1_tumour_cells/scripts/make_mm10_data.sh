wget ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz
gunzip Mus_musculus.GRCm38.93.gtf.gz
python gtf.py --attribute_name gene_name Mus_musculus.GRCm38.93.gtf mm10.data
