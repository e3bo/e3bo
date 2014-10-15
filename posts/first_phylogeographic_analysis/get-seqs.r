library('rentrez')

search <- entrez_search(db="nucleotide", term="Porcine epidemic diarrhea virus[Organism] AND (country=USA OR country=Mexico OR country=Canada)", retmax=150)
nuc <- entrez_fetch(db="nucleotide", id=search$ids, rettype="gb")
write(nuc, file='pedv-usa-mexico.gb')
system('python2 parse.py')
feats <- read.csv('pedv-usa-mexico.csv', stringsAsFactors=FALSE)

ind <- grep(':', feats$country)
feats <- feats[ind, ]
foo <- strsplit(feats$country, split=':[ ]*')
foo <- sapply(foo, '[[', 2)
foo[foo=='NorthCarolina'] <- 'North Carolina'
foo[foo=='Tennesse'] <- 'Tennessee'
foo[foo=='CO'] <- 'Colorado'
stopifnot(foo %in% state.name)
ind <- match(foo, state.name)
feats$abb <- state.abb[ind]
feats$datef <- as.Date(feats$date, '%d-%b-%Y')
test <- !is.na(feats$datef)
feats <- feats[test,]

query <- paste(feats$accession, collapse=' ')
res <- entrez_search(db='nucleotide', term=query, retmax=100)
seqs <- entrez_fetch(db='nucleotide', rettype='fasta', id=res$ids)
write(seqs, 'pedv.fasta')

library(ape)
dna <- read.dna("pedv.fasta", "fasta")

labs <- names(dna)
labs <- strsplit(labs, split='\\|')
acsn <- sapply(labs, '[[', 4)
feats$tname <- paste(feats$datef, feats$abb, feats$accession, sep='_')
ind <- match(acsn, feats$accession)
names(dna) <- feats$tname[ind]
write.dna(dna, file='renamed.fasta', format='fasta')

seqali <- muscle(read.dna("pedv.fasta", "fasta"))
