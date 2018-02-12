# digest_planner


Add plasmid sequence to plasmid.txt. Easiest way is to remove it and make a new file:

'''
rm plasmid.txt
vim plasmid.txt
'''

Then just run digestPlanner.py. It will display a list of enzymes and the expected pattern they would produce. One small caveat is that it doesn't currently wrap around the plasmid, so if there is a recognition site that starts in the last couple bases, and includes the first couple bases, it won't be found. 

Enzymes are stored as a dictionary in digest_planner.py (I will update to import from csv soon).
It only includes enzymes that are in our lab. I plan to update it to also store the buffer and whether it is an HF enzyme.

The scoring is arbitrary. It negatively scores very small/large bands, closely spaced bands, numerous bands. I plan to also include a comparison to the expected pattern for the DEST vectors.

  
