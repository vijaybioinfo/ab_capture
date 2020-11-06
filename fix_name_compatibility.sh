## Fixing prefix mess...
cd /mnt/BioAdHoc/Groups/vd-vijay/cramirez/amg319/raw/NV035
ls -loh
mkdir archive
FNAMES=(`ls -d count/*Gex`)
for FNAME in ${FNAMES[@]}; do
  echo ${FNAME}
  FIND_NAME=`basename ${FNAME} | sed -E 's/_Gex//; s|[0-9]{1,}_||'`
  CURRENT_NAME=`dirname ${FNAME}`/`ls $(dirname ${FNAME}) | grep ${FIND_NAME}_CITE`
  echo ln -s `pwd`/${CURRENT_NAME/count/archive} ${CURRENT_NAME/_Gex/_CITE}
done
