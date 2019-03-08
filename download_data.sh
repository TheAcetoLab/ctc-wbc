mkdir -p data

if [ ! -e "data/sce_hs.rds" ]; then
  curl -o data/sce_hs.rds.gz  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109761&format=file&file=GSE109761%5Fsce%5Fhs%2Erds%2Egz"
  gunzip data/sce_hs.rds.gz
fi

if [ ! -e "data/sce_mm.rds" ]; then
  curl -o data/sce_mm.rds.gz  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109761&format=file&file=GSE109761%5Fsce%5Fmm%2Erds%2Egz"
  gunzip data/sce_mm.rds.gz
fi