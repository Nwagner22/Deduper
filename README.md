# Deduper

### Usage:
To run Deduper use the command **python wagner_deduper.py** in combination with the flags below.


### Flags:
* -f: Use to specify the file path to the sam file to be de-duplicated  REQUIRED
* -o: Use to specify the extention to the input file name.  DEFAULT: _deduped.sam
* -u: Use to specify the file path to the file containing known UMIs  DEFAULT: UMI96.txt
* -p: NOT YET IMPLEMENTED
* -h: Help

### Required packages:
* re (regular expressions)
* os (operating system)
* argparse

**Example:**
```
python dedupifier.py -s Dataset3.sam 
```

#### OR

Use the **run_deduper.srun** script on cluster systems that support slurm.
```
dedupifier_path='/projects/bgmp/nwagner2/Genom_research_lab/Assignments/Deduper/Deduper'
dataset_path='/projects/bgmp/shared/deduper'

/usr/bin/time -v python $dedupifier_path/dedupifier.py -s $dataset_path/Dataset1_sorted.sam

/usr/bin/time -v python $dedupifier_path/dedupifier.py -s $dataset_path/Dataset2_sorted.sam

/usr/bin/time -v python $dedupifier_path/dedupifier.py -s $dataset_path/Dataset3_sorted.sam
```


### Runtime:

The biggest of the three datasets this was tested on was 1.6G and it took roughly 45 seconds to deduplicate it.
```
Removed 369005 duplicates
	Command being timed: "python /projects/bgmp/nwagner2/Genom_research_lab/Assignments/Deduper/Deduper/dedupifier.py -s /projects/bgmp/nwagner2/Genom_research_lab/Assignments/Deduper/Dataset1_sorted.sam -u /projects/bgmp/nwagner2/Genom_research_lab/Assignments/Deduper/Data/UMI96.txt"
	User time (seconds): 8.13
	System time (seconds): 0.28
	Percent of CPU this job got: 98%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:08.51
	Maximum resident set size (kbytes): 187752
	Exit status: 0
Removed 628602 duplicates
	Command being timed: "python /projects/bgmp/nwagner2/Genom_research_lab/Assignments/Deduper/Deduper/dedupifier.py -s /projects/bgmp/nwagner2/Genom_research_lab/Assignments/Deduper/Dataset2_sorted.sam -u /projects/bgmp/nwagner2/Genom_research_lab/Assignments/Deduper/Data/UMI96.txt"
	User time (seconds): 10.95
	System time (seconds): 0.39
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:11.43
	Maximum resident set size (kbytes): 252432
	Exit status: 0
Removed 1618545 duplicates
	Command being timed: "python /projects/bgmp/nwagner2/Genom_research_lab/Assignments/Deduper/Deduper/dedupifier.py -s /projects/bgmp/nwagner2/Genom_research_lab/Assignments/Deduper/Dataset3_sorted.sam -u /projects/bgmp/nwagner2/Genom_research_lab/Assignments/Deduper/Data/UMI96.txt"
	User time (seconds): 45.43
	System time (seconds): 1.38
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:47.21
	Maximum resident set size (kbytes): 280024
	Exit status: 0
```