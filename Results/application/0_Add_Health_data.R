# Data
rm(list = ls())

library(readstata13)
library(ggplot2)
library(CDatanet)
library(DataCombine)
library(stargazer)
library(dplyr)

filname        <- "MydataCount.rda"
if (file.exists(filname)) {
  load(file = filname)
} else {
  mydata       <- read.dta13("cleandta.dta")
  mydata       <- mydata[order(mydata$sschlcde),]
  
  mislist      <- c(55555555,77777777,88888888,99999999,99959995)
  f_coln       <- c(paste0("mf", 1:5, "aid"), paste0("ff", 1:5, "aid"))
  
  all.va.names <- c(
    "age",        "male",        "female",      "grade",       "hispanic",
    "racewhite",  "raceblack",   "raceasian",   "raceother",   "yearinschl",  "withmom",
    "melhigh",    "mehigh",      "memhigh",     "memiss",      "mjhome",      "mjprof",
    "mjother",    "mjmiss",    "withdad",       "withbothpar", "felhigh",     "fehigh",
    "femhigh",    "femiss",      "fjhome",      "fjprof",      "fjother",     "fjmiss",
    "nclubs") 
  
  va.names      <- c("nclubs", "age", "male", "hispanic", "raceblack",
                     "raceasian", "raceother", "yearinschl", "withbothpar", "melhigh", "memhigh",
                     "memiss", "mjprof", "mjother",  "mjmiss")
  
  apply(mydata[, va.names], 2, function(w) sum(is.na(w)))
  
  
  keep1         <- as.logical((rowSums(is.na(mydata[,va.names])) == 0))
  mydata        <- mydata[keep1,]
  
  
  
  # remove friend from different groups
  # remove self friendship
  # remove friend non found
  N             <- nrow(mydata)
  repf          <- rep(0,N)
  dscf          <- rep(0,N)
  sfre          <- rep(0,N)
  nffr          <- rep(0,N)
  for (i in 1:N) {
    jcont       <- 0
    for (j in f_coln) {
      k         <- which(mydata$aid == mydata[i, j])
      
      if (length(k) != 0) {
        # remove if already listed
        if (jcont > 0) {
          if(mydata[i, j] %in% mydata[i, f_coln][1:jcont]) {
            mydata[i, j]   <- 0
            repf[i]        <- repf[i] + 1
          }
        }
        # remove if different school
        if(mydata[i, "sschlcde"] != mydata[k, "sschlcde"]) {
          mydata[i, j]   <- -1
          dscf[i]        <- dscf[i] + 1
        }
        # remove if self frienship
        if(mydata[i, "aid"] == mydata[k, "aid"]) {
          mydata[i, j]   <- -2
          sfre[i]        <- sfre[i] + 1
        }
      }
      else {
        if (!((mydata[i, j] %in% mislist) | is.na(mydata[i, j]))) {
          mydata[i, j]   <- -3
          nffr[i]        <- nffr[i] + 1
        }
      }
      jcont       <- jcont + 1
    }
  }
  
  cat("remove", sum(repf), "link(s) because links repeated: their code are recoded as 0", "\n")
  cat("remove", sum(dscf), "link(s) from different schools: their code are recoded as -1", "\n")
  cat("remove", sum(sfre), "self-friendship(s): their code are recoded as -2", "\n")
  cat("remove", sum(nffr), "non-found friends: their code are recoded as -3", "\n")
  rm(list = c("i", "j", "k"))
  
  # Compute number of friends
  mf_coln   <- paste0("mf", 1:5, "aid")  
  ff_coln   <- paste0("ff", 1:5, "aid")  
  nfriends  <- apply(mydata[,f_coln], 1, function(w) sum(!(w %in% c(0, -1, -2, -3, NA, mislist)))) 
  nmfriends <- apply(mydata[,mf_coln], 1, function(w) sum(!(w %in% c(0, -1, -2, -3, NA, mislist)))) 
  nffriends <- apply(mydata[,ff_coln], 1, function(w) sum(!(w %in% c(0, -1, -2, -3, NA, mislist))))
  
  # remove schools wihout link or schools that lose more link than they have
  tmp       <- cbind(aggregate(data.frame("friend"   = nfriends,
                                    "diffsch"  = dscf,
                                    "self"     = sfre,
                                    "nofound"  = nffr), 
                         list("IDsch" = mydata$sschlcde), 
                         FUN = "mean"),
                     "N" = as.numeric(table(mydata$sschlcde)))
  scidlist  <- tmp[(tmp$friend > tmp$nofound) & tmp$friend > 0 & tmp$N > 100, "IDsch"]
  keep2     <- mydata$sschlcde %in% scidlist
  
  
  mydata    <- mydata[keep2,]
  dscf      <- dscf[keep2]
  sfre      <- sfre[keep2]
  nffr      <- nffr[keep2]
  repf      <- repf[keep2]
  
  #remove schools with less than 100
  
  N         <- nrow(mydata)
  
  # y
  y         <- mydata$nclubs
  
  # X
  X         <- as.matrix(mydata[,va.names[-1]])
  
  
  # friends
  aid       <- mydata$aid
  
  igroup1   <- aggregate(data.frame("min" = 0:(N-1)), 
                         list("IDsch" = mydata$sschlcde), 
                         FUN = "min")
  
  igroup2   <- aggregate(data.frame("max" = 0:(N-1)), 
                         list("IDsch" = mydata$sschlcde), 
                         FUN = "max")
  
  igroup    <- cbind(igroup1$min, igroup2$max)
  ngroup    <- nrow(igroup)
  
  Gnet      <- list()
  expv      <- c("age", "male", "hispanic", "raceblack", "raceasian", "raceother", "yearinschl", "withbothpar", "melhigh", "memhigh",
                     "memiss", "mjprof", "mjother",  "mjmiss")
  va.net    <- c("age", "male", "hispanic", "racewhite", "raceblack", "raceasian", "yearinschl", "mjprof")
  nvar      <- length(va.net)
  
  typevar   <- rep(3, nvar)
  typevar[va.net %in% c("age", "yearinschl")]  <- 1
  typevar[va.net %in% "male"]                  <- 2
  
  datanet   <- matrix(0, 0, nvar)
  
  for (k in 1:ngroup) {
    # Network
    cat("group ", k, "/", ngroup, "\n")
    n1      <- igroup[k,1] + 1
    n2      <- igroup[k,2] + 1
    tmpfr   <- apply(mydata[n1:n2, f_coln], 1, function(w) which(aid[n1:n2] %in% w))
    Nk      <- n2 - n1 + 1
    
    Gk      <- matrix(data = 0L, nrow = Nk, ncol = Nk)
    
    for (i in 1:Nk) {
      Gk[i, tmpfr[[i]]] <- 1  
    }
    
    Gnet[[k]]  <- Gk
    
    # Network formation
    Xnet           <- list()
    for (j in 1:nvar) {
      Xnetj        <- mydata[n1:n2, va.net[j]]
      
      if (typevar[j] == 1) {
        Xnetj      <- kronecker(Xnetj, t(Xnetj), function(x, y) abs(x - y))
      } else {
        if (typevar[j] == 2) {
          Xnetj    <- kronecker(Xnetj, t(Xnetj), function(x, y) 1*(x == y))
        } else {
          Xnetj    <- kronecker(Xnetj, t(Xnetj), function(x, y) 1*((x == 1) & (y == 1)))
        }
      }
      
      diag(Xnetj)  <- NA
      Xnetj        <- c(Xnetj)
      Xnet[[j]]    <- Xnetj[!is.na(Xnetj)]
    }
    Xnet           <- do.call("cbind", Xnet)
    datanet        <- rbind(datanet, Xnet)
  }
  
  # Compute number of misslink
  miss1  <- apply(mydata[,f_coln], 1, function(w) sum((w == mislist[1]), na.rm = TRUE)) 
  miss2  <- apply(mydata[,f_coln], 1, function(w) sum((w == mislist[2]), na.rm = TRUE)) 
  miss3  <- apply(mydata[,f_coln], 1, function(w) sum((w == mislist[3]), na.rm = TRUE))  
  miss4  <- apply(mydata[,f_coln], 1, function(w) sum((w == mislist[4]), na.rm = TRUE))  
  miss5  <- apply(mydata[,f_coln], 1, function(w) sum((w == mislist[5]), na.rm = TRUE))  
  
  mydata        <- cbind(mydata, 
                         "nfriend"     = nfriends[keep2],
                         "nffriend"    = nffriends[keep2],
                         "nmfriend"    = nmfriends[keep2],
                         "mis55555555" = miss1,
                         "mis77777777" = miss2,
                         "mis88888888" = miss3,
                         "mis99999999" = miss4,
                         "mis99959995" = miss5)
  
  # schools
  school        <- aggregate(x   = data.frame(size = mydata$aid),
                             by  = list("IDsch"    = mydata$sschlcde),
                             FUN = length)
  
  tmp           <- aggregate(x   = data.frame(friends  = mydata$nfriend,
                                              mfriends = mydata$nmfriend,
                                              ffriend  = mydata$nffriend,
                                              nclubs   = mydata$nclubs),
                             by  = list("IDsch" = mydata$sschlcde),
                             FUN = mean)
  
  school        <- cbind(school, tmp[,-1])
  school        <- school[order(school$size),]
  n.school      <- nrow(school)
  school$school <- 1:n.school
  
  
  Xd            <- matrix(0, N, n.school); 
  colnames(Xd)  <- paste0("I_", 1:n.school)
  for (i in 1:N) {
    Xd[i,which(scidlist == mydata$sschlcde[i])] <- 1
  }
  mydata            <- cbind(mydata, Xd)
  colnames(datanet) <- va.net
  datanet           <- as.data.frame(datanet)
  
  rm(list = ls()[!(ls() %in% c("Gnet", "mydata", "school", "n.school", "repf",
                               "dscf", "sfre", "nffr", "filname", "expv", "datanet", "va.net"))])
  
  eff           <- paste0("I_", 1:n.school)
  GX            <- peer.avg(norm.network(Gnet), mydata[,expv]); colnames(GX) <- paste0("g", expv)
  mydata        <- cbind(mydata, GX)
  expvc         <- c(expv, paste0("g", expv))
  save(list = ls(all = TRUE), file = filname)
}

# Descriptive statistics
my.sum <- function(x) {
  out <- c(mean(x, na.rm = TRUE),
           sd(x, na.rm = TRUE),
           min(x, na.rm = TRUE),
           quantile(x, 0.25, na.rm = TRUE),
           median(x, na.rm = TRUE),
           quantile(x, 0.75, na.rm = TRUE),
           max(x, na.rm = TRUE))
  names(out) <- c("Mean", "Sd.", "Min", "1st Qu.", "Median", "3rd Qu.", "Max")
  return(out)
}

mydata       <- mydata %>% mutate(nclub_male = ifelse(male == 1, nclubs, NA),
                                  nclub_female = ifelse(male == 0, nclubs, NA))

all.va.names <- c( "age", "male", "hispanic", "racewhite", "raceblack", "raceasian", "raceother",   
                   "yearinschl", "withbothpar", "melhigh", "mehigh", "memhigh", "memiss", "mjhome",      
                   "mjprof", "mjother", "mjmiss", "nclubs", "nclub_male", "nclub_female") 

allvar <- mydata[,all.va.names]

cat("Data Summary", "\n")
desc.data <- t(apply(allvar, 2, my.sum))


c1.names  <- c("Age", "Male", "Hispanic", "Race", rep("", 4), "Year in school",
               "With both parents", "Mother Education", rep("", 4), "Mother Job", 
               rep("", 4), "Number of activities", rep("", 2))
c2.names  <- c(rep("", 4), "White", "Black", "Asian", "Other", rep("", 3), "high",
               "< high", "> high", "Missing", "", "Stay home",  "Professional", "Other",
               "Missing", "", "For males", "For females")

desc.data <- as.data.frame(cbind(format(round(desc.data[,c(1,2)],3), 3), desc.data[,-c(1,2)]))

for (i in 1:7) {
  desc.data[,i] <- as.character(desc.data[,i])
}

desc.data <- InsertRow(desc.data, NewRow = rep("NA", 7), RowNum = 4)
desc.data <- InsertRow(desc.data, NewRow = rep("NA", 7), RowNum = 11)
desc.data <- InsertRow(desc.data, NewRow = rep("NA", 7), RowNum = 16)

desc.data <- as.data.frame(cbind(c1.names, c2.names, desc.data))   

for (i in 3:9) {
  desc.data[grepl("NA", desc.data[,i]),i]     <- ""
}
rownames(desc.data) <- NULL
# Table D.2: Data Summary
print(desc.data)
stargazer(desc.data, summary = FALSE, rownames = FALSE)

# Figure 3: Distribution of the number of extracurricular activities
(graph  <- ggplot(data = data.frame(y = mydata$nclubs), aes(x = y)) +
    geom_bar(color = "black", fill = "#eeeeee") + 
    theme_bw() + xlab("") + ylab("Frequency")  + 
    theme(strip.text = element_text(face = "italic"), 
          text = element_text(size = 12, family = "Palatino"),
          axis.title = element_text(size = 12, family = "Palatino")))

ggsave("Nactivities.pdf", path = "_output", plot = graph, device = "pdf", width = 6, height = 2.5)