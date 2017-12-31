# DATA=lastfm
DATA=dblp
# DATA=diggs
# DATA=flixster
ITEMSET_TRAIN='../data/'$DATA'.train.itemset'
source ~/pyenv/taim/bin/activate

for p in 1 5 10 15 20
do
    UITEMSAMPLE='../data/'$DATA'.'$p'p.up.sample'
    LITEMSAMPLE='../data/'$DATA'.'$p'p.low.sample'
    N=`expr $p \* 80`
    echo "python kmeans.py $ITEMSET_TRAIN $LITEMSAMPLE $UITEMSAMPLE -i 100 -n $N -m 10"
    python kmeans.py $ITEMSET_TRAIN $LITEMSAMPLE $UITEMSAMPLE -i 100 -n $N -m 10
done
