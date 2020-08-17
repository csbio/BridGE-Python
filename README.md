# BridGE-Python

This is the github repository for the Python version of BridGE


## Data Processing
To run:
`python3 bridge.py --job=DataProcess --options`

Options:
--plinkFile = path to plink input file without extention. 
Example:
`python3 bridge.py --job=DataProcess --plinkFile=data/gwas_data_final`

## Compute Interaction network
To run:
`python3 bridge.py --job=ComputeInteraction --options`

Options:
--model = model, RR - RD - DD - combined
--nWorker = number of parallel workers
--samplePerms = number of sample permuations
Example:
`python3 bridge.py --job=ComputeInteraction --model=RR --nWorker=6 --samplePerms=10`

## Sample Permutation
To run:
`python3 bridge.py --job=ComputeInteraction --options`

Options:
--model = model, RR - RD - DD - combined
--nWorker = number of parallel workers
--samplePerms = number of sample permuations
--snpPerms = number of snp permutations used in generatin pathway stats
--binaryNetwork = Binary flag to use binarized interactions
