# digest_planner

You may need to install tabulate for the script to work.

```
pip install tabulate
```

Add plasmid sequence to plasmid.txt. Easiest way is to remove it and make a new file:

```
rm plasmid.txt
vim plasmid.txt
```

Then just run digestPlanner.py. It will display a list of enzymes and the expected pattern they would produce.

Enzymes are stored as a csv file. It only includes enzymes that are in our lab, but may not be completely up to dat. I plan to update it to also store the buffer and whether it is an HF enzyme.

The scoring is arbitrary. It penalizes very small/large bands, closely spaced bands, and numerous bands.

Two small caveats:

1) The plasmid is treated as a linear sequence, not a circular sequence. Therefore if there is a recognition sequence that includes bases before and after the first base of the seqeunce, it will not be found. I don't expect this issue to be common. If you run into it, re-index your plasmid on Benchling. I can implement a more general fix if needed.

2) The cutting locations for the enzyme aren't saved. The band pattern is based on the location of the recognition sequence, not where the cut actually occurs. The band patterns are rounded to the nearest hundred and it's possible that this error will cause a bands displayed size to be as much as 100bp larger or smaller than the expected. I think the effort required to fix this is not worth it.
