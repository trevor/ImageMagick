SRCDIR=`dirname $0`
SRCDIR=`cd $SRCDIR && pwd`
TOPSRCDIR=`cd $srcdir && pwd`
cd tests || exit 1
PASSPHRASE="${TOPSRCDIR}/utilities/tests/passphrase.txt"
