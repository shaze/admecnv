.mode tabs
CREATE TABLE delly_cnv(chrom text, start int, end int, cn0, cn1 int, cn2 int, cn3 int, cn4 int, cn5
int, cn6 int, cn7plus int);
create table variants ( chrom  text, start  int, end int, type  text, ac  int, ah  text, classifier
text);
.import manta.tsv variants
.import delly.tsv variants
.import smoove.tsv variants
.import gridss.tsv variants
.import "delly_cnv.tsv" delly_cnv
create index dcnv_end on delly_cnv(chrom,end);
create index dcnv_start on delly_cnv(chrom,start);
create index classifier_ind on variants(classifier);
create index start_class_ind  on variants(chrom, start, classifier);
create index start_class_ind  on variants(chrom, end, classifier)
create index start_class_ind  on variants(chrom, start, end, classifier)


.quit
