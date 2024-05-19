# Summary results

Currently this contains summary results from all the African data we have used in the paper.

We will work with the data providers to provide regional summary data, consistent with ethics.

## Calls from individual tools

### Structural variants

For
  delly.tsv
  manta.tsv
  gridss.tsv
  smoove.tsv


The TSV files provide summaries for the calss. The fields are
-- chromosome
-- start and end position
-- type of variant
-- ac : number of calls made
-- ah : number of homozygous calls made

As an example, if there were 10 individuals who had homozygous calls and 100 who had heterozygous calls, then ac will be 120 and ah will be 20

The exception to this is gridss where ac is the number of individuals for whom an SV was called


### Copy number variants

`delly_cnv.tsv` shows the results from delly in CNV mode

The fields are
chrom -- chromosome
start
end
cn0  -- number of calls with 0 copies
cn1  -- number of calls with 1 copy
cn2  -- number of calls with 2 copies (expected in a diploid organism)
cn3  -- number of calls with 3 copies
cn4  -- number of calls with 4 copies
cn5  -- number of calls with 5 copies
cn6  -- number of calls with 6 copies
cn7plus -- number of calls with 7 or more copies