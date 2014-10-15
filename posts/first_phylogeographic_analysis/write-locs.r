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

locs <- data.frame(abb=unique(feats$abb))
ind <- match(locs$abb, state.abb)
locs$lat <- state.center$y[ind]
locs$long <- state.center$x[ind]
colnames(locs) <- NULL
write.table(locs, file='locations.csv', sep='\t', row.names=FALSE, quote=FALSE)

