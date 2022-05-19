#!/bin/bash

# odgi_ sort -i DRB1-3123_unsorted.og --threads 2 -P -Y -o DRB1-3123_sorted.og  -u DRB1-3123_sorted.og_snapshot
# for OG in *_snapshot*; do odgi viz -i "$OG" -o "$OG".png; done
# this gives us the highest height in pixels of all our candidates
# identify *snapshot*.png | cut -f 3 -d " " | cut -f 2 -d "x" | sort | tail -n 1
# identify -format '%h\n' *snapshot*.png | sort | tail -n 1

# convert -size 1744x571 xc:white DRB1-3123_sorted.og_snapshot29.png -gravity North -composite test.png

#HEIGHT=`identify -format '%h\n' *snapshot*.png | sort -n | tail -n 1`

#rm -rf *snapshot*.resized.png

#for PNG in *snapshot*.png
#do
    # https://legacy.imagemagick.org/discourse-server/viewtopic.php?t=11525
#    convert -size 1744x"$HEIGHT" xc:white "$PNG" -gravity North -composite "$PNG".resized.png
#done

# https://askubuntu.com/questions/648244/how-do-i-create-an-animated-gif-from-still-images-preferably-with-the-command-l
#convert -delay 0 -dispose previous -loop 0 `ls -v *snapshot*.resized.png` DRB1-3123_sorted.png DRB1-3123_sorted_snapshots.gif

# odgi layout -i DRB1-3123_unsorted.og -o DRB1-3123_unsorted.og.lay -P --threads 2 -u DRB1-3123_unsorted.og.lay.snapshot
# for OG in *lay.snapshot*; do odgi draw -i DRB1-3123_unsorted.og -c "$OG" -p "$OG".png -H 500 -C -w 10; done

#WIDTH=`identify -format '%w\n' *lay.snapshot*.png | sort -n | tail -n 1`

#rm -rf *lay.snapshot*.resized.png

#for PNG in *lay.snapshot*.png
#do
    # https://legacy.imagemagick.org/discourse-server/viewtopic.php?t=11525
#    convert -size "$WIDTH"x500 xc:white "$PNG" -gravity Center -composite "$PNG".resized.png
#done

#convert -delay 0 -dispose previous -loop 0 `ls -v *lay.snapshot*.resized.png` DRB1-3123_sorted.lay_snapshots.gif
