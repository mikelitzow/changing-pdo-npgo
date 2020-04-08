# download PDO/NPGO and process into winter means

# load pdo
download.file("http://jisao.washington.edu/pdo/PDO.latest", "data/pdo") 
names <- read.table("data/pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("data/pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>% # not loading tidyr because I think it conflicts with maps!
  arrange(YEAR)

pdo$win.yr <- ifelse(pdo$month %in% c("NOV", "DEC"), pdo$YEAR+1, pdo$YEAR)

pdo <- pdo %>%
  filter(month %in% c("NOV", "DEC", "JAN", "FEB", "MAR"))

pdo <- tapply(pdo$value, pdo$win.yr, mean)

pdo <- data.frame(year=names(pdo), PDO=pdo)

pdo <- pdo %>%
  filter(year %in% 1951:2018)

write.csv(pdo, "SOMscripts/InputData/PDO.csv", row.names = F)

# load npgo
download.file("http://www.oces.us/npgo/npgo.php", "data/npgo")
npgo <- read.table("data/npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))

npgo$win.yr <- ifelse(npgo$month %in% 11:12, npgo$Year+1, npgo$Year)

npgo <- npgo %>%
  filter(month %in% c(11,12,1:3))

npgo <- tapply(npgo$value, npgo$win.yr, mean)

npgo <- data.frame(year=names(npgo), NPGO=npgo)

npgo <- npgo %>%
  filter(year %in% 1951:2018)

write.csv(npgo, "SOMscripts/InputData/NPGO.csv", row.names = F)
