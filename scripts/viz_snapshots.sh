# Example usage:
# ./viz_snapshots "DRB1-3123.fa.gz.gfa.og.Y.og_*"
FILES=$1
for f in $FILES; do odgi viz -i ${f} -o ${f}.png; done;