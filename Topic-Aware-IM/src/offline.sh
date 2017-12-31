DATA=lastfm
TEDGES='../data/'$DATA'.tedge'
LIST='../data/'$DATA'.list'

bound1=100
bound2=10  # out neighbor bound
k=50
UITEMSAMPLE='../data/'$DATA'.up.sample'
LITEMSAMPLE='../data/'$DATA'.low.sample'
UBOUND='../data/'$DATA'.up.bound'
LBOUND='../data/'$DATA'.low.bound'

make clean && make offline
./offline -l $TEDGES $bound1 $bound2 $LIST
# ./offline -us $TEDGES $LIST $bound1 $k $UITEMSAMPLE $UBOUND
# ./offline -ls $TEDGES $LIST $bound1 $k $LITEMSAMPLE $LBOUND
