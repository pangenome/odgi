import argparse

from kneed import KneeLocator
from math import floor

import pyranges as pr
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = [15, 10]


parser = argparse.ArgumentParser()

parser.add_argument(
    "-i",
    "--bed-in",
    dest="bed_in",
    help="Input path to BED file.",
    required=True
)

parser.add_argument(
    "-S",
    "--kneed-S",
    dest="kneed_S",
    help="Kneed S parameter https://kneed.readthedocs.io/en/stable/parameters.html#s.",
    required=True
)

parser.add_argument(
    "-o",
    "--prefix-out",
    dest="prefix_out",
    help="Output prefix path for BED and PNG file(s).",
    required=False
)

args = parser.parse_args()

INPUT=args.bed_in
S=args.kneed_S
if args.prefix_out is not None:
    OUTPUT_PREFIX=args.prefix_out
else:
    OUTPUT_PREFIX=INPUT

# debugging
#print(INPUT)
#print(S)
#print(OUTPUT_PREFIX)

def truncate(f, n):
    return floor(f * 10 ** n) / 10 ** n

# TODO What if we have several S?
bed = pr.read_bed(INPUT, as_df=True)
y = bed["Strand"]
x = range(1, len(y) + 1)

for kneed_s in str(S).split(","):

    kn = KneeLocator(x, y, curve='convex', direction='decreasing', S=float(kneed_s))
    # print(kn.knee)

    plt.xlabel('bin')
    plt.ylabel('ratio')
    plt.title(INPUT.split('/')[-1] + ': ' + ' S=' + str(kneed_s) + ', knee=' + str(kn.knee) + ', min_ratio=' + str(round(y[kn.knee-1], 2)))
    plt.plot(x, y, 'bo-')
    plt.vlines(kn.knee, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')

    PNG_OUT = OUTPUT_PREFIX + "_S" + kneed_s + ".png"
    BED_OUT = OUTPUT_PREFIX + "_S" + kneed_s + ".bed"

    plt.savefig(PNG_OUT)

    bed_kneed = bed[:kn.knee]
    bed_kneed.to_csv(BED_OUT, sep="\t", index=False, header=False)
    plt.close()
    plt.cla()
    plt.clf()