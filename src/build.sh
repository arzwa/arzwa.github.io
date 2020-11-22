if [ -z "$2" ]; then
    OUT=`basename $1 .md`.html
    echo "No <out> argument supplied, taking $OUT"
else
    OUT=$2
fi

pandoc -o $OUT $1 --mathml --css ./css/bb.css --standalone --metadata title=" "
