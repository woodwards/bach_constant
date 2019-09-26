# read data files into data frames

# load run control file (put it in run dir so it gets saved with output)
file_name <- paste(out_path, 'bachrunlist.dat', sep='')
runlist <- read_tsv(file=file_name, skip=0, col_types=cols(),
                    col_names=c('runname', 'catchfile', 'optfile', 'control', 'nits',
                                'chem1i', 'chem1wt', 'chem1ae', 'chem1re',
                                'chem2i', 'chem2wt', 'chem2ae', 'chem2re'))

# add catchment names and ordering information
runlist <- runlist %>%
  mutate(
    catchname=rep(c('Tahunatara', 'Otamakokore', 'Waiotapu', 'Pokaiwhenua', 'Piako', 'Waitoa', 'Waihou', 'Puniu'), times=3),
    shortname=rep(c('Ta', 'Ot', 'Wp', 'Po', 'Pi', 'Wt', 'Wh', 'Pu'), times=3),
    catchseq=rep(c(5, 4, 2, 3, 8, 7, 1, 6), times=3),
    period=rep(c('a', 'b', 'c'), each=8),
    periodseq=rep(c(1, 2, 3), each=8),
    setname=paste(shortname, '(', period, ')', sep=''),
    setseq=periodseq+(catchseq-1)*3
    )
    
# load options files
options <- tibble() 
dir_path <- 'data/'
files <- tibble(fname=list.files(path=dir_path, pattern='\\options..dat$'), path=dir_path)
for (i in 1:nrow(files)){
  file_name <- paste(files$path[i], files$fname[i], sep='')
  temp <- read_tsv(file_name, skip=0, col_names=c(files$fname[i]), col_types=cols())
  options <- rbind(options, t(temp))
}
names(options) <- c('nchem', 'startrun', 'startcalib', 'endcalib', 'startvalid', 'endvalid')

# load observed data files
data <- tibble() 
dir_path <- 'data/'
files <- tibble(fname=list.files(path=dir_path, pattern='\\_data.dat$'), path=dir_path)
for (i in 1:nrow(files)){
  file_name <- paste(files$path[i], files$fname[i], sep='')
  temp <- read_tsv(file_name, 
                   na="NaN", skip=0, progress=FALSE, n_max=Inf, col_types=cols(),
                   col_names=c('date', 'flow', 
                               'Ca', 'Mg', 'K', 'Na', 'Zn', 'HCO3', 'SO4', 'Cl',
                               'Si', 'EC', 'DRP', 'NDP', 'TP', 'NNN', 'NH4', 'TN'))
  keeps <- c('date', 'flow', 'TP', 'TN') # columns to keep
  temp <- temp[keeps] %>%  
    mutate(
      file=files$fname[i], # add file name
      area=c(45, 103, 432, 519, 208, 802, 298, 121)[i],
      shortname=c('Ot', 'Pi', 'Po', 'Pu', 'Ta', 'Wh', 'Wp', 'Wt')[i],
      catchseq=c(4, 8, 3, 6, 5, 1, 2, 7)[i],
      flow2=flow/area*86.4, # mm
      period=ifelse(date>=37257 & date<=39082,'02-06',
                    ifelse(date>=39083 & date<=40908,'07-11',
                           ifelse(date>=40909 & date<=42735,'12-16',
                                  'None'))),
      period2=ifelse(date>=37257 & date<=39082,'a',
                     ifelse(date>=39083 & date<=40908,'b',
                            ifelse(date>=40909 & date<=42735,'c',
                                   'None'))),
      periodseq=ifelse(date>=37257 & date<=39082,1,
                     ifelse(date>=39083 & date<=40908,2,
                            ifelse(date>=40909 & date<=42735,3,
                                   0))),
      setname=paste(shortname, '(', period2, ')', sep=''),
      rdate=as.Date(date, origin="1899-12-30"),
      year=lubridate::year(rdate),
      setseq=periodseq+(catchseq-1)*3
      ) 
  data <- rbind(data, temp)
}
