setwd("C:/Users/Sagar/OneDrive/Documents/EI/HW4")
library(lpSolveAPI)
library(ggplot2)
library(knitr)
library(rmarkdown)

#Load data
demand.mw <- read.csv("demand.csv",stringsAsFactors = FALSE) 
gen.info <- read.csv("gen.info.csv") 
import.limits <- read.csv("import.limits.csv") 
wind.sites <- read.csv("wind.sites.csv") 
sites.wind.power.per <- read.csv("sites.wind.power.per.csv",stringsAsFactors = FALSE) 
trans.limits <- read.csv("trans.limits.csv") 
rated.cap.per.mw <- 3

#I define a few objects that will make things easier to generalize later
n.sites <- nrow(wind.sites)
n.hrs <- nrow(demand.mw)
n.regions <- nrow(gen.info)

#Set transmission limits
t.west_east.limit.mw <- trans.limits$trans.limit.mw[which(trans.limits$from == 
                                                            "west" & trans.limits$to == "east")]
t.east_downstate.limit.mw <- trans.limits$trans.limit.mw[which(trans.limits$from == 
                                                                 "east" & trans.limits$to == "downstate")]

#Calculate baseload generation
gen.base.west.mw <- sum(gen.info[which(gen.info$region == "west"), c("ann.gen.coal.gwh", 
                                                                     "ann.gen.nuc.gwh", "ann.gen.hydro.gwh")]) * 1000/n.hrs
gen.base.east.mw <- sum(gen.info[which(gen.info$region == "east"), c("ann.gen.coal.gwh", 
                                                                     "ann.gen.nuc.gwh", "ann.gen.hydro.gwh")]) * 1000/n.hrs
gen.base.downstate.mw <- sum(gen.info[which(gen.info$region == "downstate"), 
                                      c("ann.gen.coal.gwh", "ann.gen.nuc.gwh", "ann.gen.hydro.gwh")]) * 1000/n.hrs
gen.base.total.mw <- sum(gen.base.west.mw, gen.base.east.mw, gen.base.downstate.mw)
data.frame(gen.base.west.mw, gen.base.east.mw, gen.base.downstate.mw, gen.base.total.mw)

#Find limits for 'other' generation
cap.gen.other.west.mw <- gen.info$cap.gen.nonbase.mw[which(gen.info$region == 
                                                            "west")] + sum(import.limits[which(import.limits$region == "west"), -1])
cap.gen.other.east.mw <- gen.info$cap.gen.nonbase.mw[which(gen.info$region == 
                                                             "east")] + sum(import.limits[which(import.limits$region == "east"), -1])
cap.gen.other.downstate.mw <- gen.info$cap.gen.nonbase.mw[which(gen.info$region == 
                                                                  "downstate")] + sum(import.limits[which(import.limits$region == "downstate"), 
                                                                                      -1])
cap.gen.other.total.mw <- sum(cap.gen.other.west.mw, cap.gen.other.east.mw, 
                              cap.gen.other.downstate.mw)
data.frame(cap.gen.other.west.mw, cap.gen.other.east.mw, cap.gen.other.downstate.mw, 
           cap.gen.other.total.mw)

#Store sites
sites.west <- wind.sites$site[which(wind.sites$region == "west")]
sites.east <- wind.sites$site[which(wind.sites$region == "east")]
sites.downstate <- wind.sites$site[which(wind.sites$region == "downstate")]
sites <- wind.sites$site
n.sites.west <- length(sites.west)
n.sites.east <- length(sites.east)
n.sites.downstate <- length(sites.downstate)

#Store installation capacities
capacities.mw <- 15000
results.cf <- as.data.frame(capacities.mw)
names(results.cf) <- "wind.capacity.mw"

#Create array to store turbine distribution
results.turbine.distr <- as.data.frame(wind.sites)

#Find total demand
demand.mw$total.mw <- rowSums(demand.mw[, c("west", "east", "downstate")])

#CASE 1

#Set transmission limits
t.west_east.limit.mw <- trans.limits$trans.limit.mw[which(trans.limits$from == "west" & trans.limits$to == "east")]
t.east_downstate.limit.mw <- trans.limits$trans.limit.mw[which(trans.limits$from == "east" & trans.limits$to == "downstate")]

#Idenitify the number of variables
n.variables <- n.sites + n.regions*n.hrs + (n.regions-1)*n.hrs

#Number of constraints
n.cons <- 1 + n.regions*n.hrs

#Number of bounded variables 
n.bounds <- n.sites + n.regions*n.hrs + (n.regions-1)*n.hrs

#Create a for loop to evaluate different levels of wind penetration
for (j in 1:length(capacities.mw)) {
  wind.cap.total.mw <- capacities.mw[j]
  
  #Define total number of turbines
  N.t <- wind.cap.total.mw/rated.cap.per.mw
  
  #Develop MILP
  
  #Create linear program model with number of rows=num of constraints and number of columns = number of variables
  lp.case1 <- make.lp(nrow = n.cons, ncol = n.variables)
  
  #Specify Min
  lp.control(lp.case1, sense = "min")
  
  #Define coefficients for the objective function
  set.objfn(lp.case1, obj = c(rep(0, n.sites), rep(1, n.regions * n.hrs), 
                              rep(0, (n.regions - 1) * n.hrs)))
  
  #Force the number of turbines at each site to be an integer
  set.type(lp.case1, c(1:n.sites), type = "integer")
  
  #Add constraints by row 
  row.add.mode(lp.case1, "on")
  
  #Sum of number of turbines at all sites equals total number of turbines
  add.constraint(lp.case1, xt = rep(1, n.sites), type = "=", rhs = N.t, indices = c(1:n.sites))
  
  #Energy balance for West region
  for (i in 1:n.hrs) {
    add.constraint(lp.case1, xt = c(sites.wind.power.per[i, 1 + sites.west], 
                                    1, -1), type = ">=", rhs = demand.mw[i, "west"] - gen.base.west.mw, 
                   indices = c(sites.west, (n.sites + i), (n.sites + n.regions * n.hrs + 
                                                             i)))
  }
  
  #Energy balance for East region
  for (i in 1:n.hrs) {
    add.constraint(lp.case1, xt = c(sites.wind.power.per[i, 1 + sites.east], 
                                    1, 1, -1), type = ">=", rhs = demand.mw[i, "east"] - gen.base.east.mw, 
                   indices = c(sites.east, (n.sites + n.hrs + i), (n.sites + n.regions * 
                                                                     n.hrs + i), (n.sites + n.regions * n.hrs + n.hrs + i)))
  }
  
  #Energy balance for Downstate region
  for (i in 1:n.hrs) {
    add.constraint(lp.case1, xt = c(sites.wind.power.per[i, 1 + sites.downstate], 
                                    1, 1), type = ">=", rhs = demand.mw[i, "downstate"] - gen.base.downstate.mw, 
                   indices = c(sites.downstate, (n.sites + 2 * n.hrs + i), (n.sites + 
                                                                              n.regions * n.hrs + n.hrs + i)))
  }
  
  #Turn off row entry mode
  row.add.mode(lp.case1, "off")
  
  #Set bounds on decision variables
  set.bounds(lp.case1, lower = rep(0, n.variables), upper = c(wind.sites$n.t.max, 
                                                              rep(cap.gen.other.west.mw, n.hrs), rep(cap.gen.other.east.mw, n.hrs), 
                                                              rep(cap.gen.other.downstate.mw, n.hrs), rep(t.west_east.limit.mw, n.hrs), 
                                                              rep(t.east_downstate.limit.mw, n.hrs)))
  
  #Solve optimization problem
  solve(lp.case1)
  turbine.distr <- get.variables(lp.case1)[1:n.sites]
  othergen.total.mwh <- get.objective(lp.case1)
  t.west_east.total.mw <- sum(get.variables(lp.case1)[(n.sites + n.regions * 
                                                         n.hrs + 1):(n.sites + n.regions * n.hrs + n.hrs)])
  t.east_downstate.total.mw <- sum(get.variables(lp.case1)[(n.sites + n.regions * 
                                                              n.hrs + n.hrs + 1):(n.sites + n.regions * n.hrs + 2 * n.hrs)])
  cf.wind <- (sum(demand.mw$total.mw) - othergen.total.mwh - gen.base.total.mw * 
                n.hrs)/(wind.cap.total.mw * 8760)
  
  #Save the turbine distribution, othergen.total and transmission totals
  assign(sprintf("turbine.distr.case1.%s", wind.cap.total.mw), turbine.distr)
  assign(sprintf("othergen.total.mwh.case1.%s", wind.cap.total.mw), othergen.total.mwh)
  assign(sprintf("t.west_east.total.mw.case1.%s", wind.cap.total.mw), t.west_east.total.mw)
  assign(sprintf("t.east_downstate.total.mw.case1.%s", wind.cap.total.mw), 
         t.east_downstate.total.mw)
  assign(sprintf("cf.wind.case1.%s", wind.cap.total.mw), cf.wind)
}

#Add capacity factors to results.cf data frame
results.cf$case1 <- (cf.wind.case1.15000)

#Add wind turbine distribution 
results.turbine.distr$n.t.case1.15gw <- c(turbine.distr.case1.15000)

#CASE 2

#Create new objects related to cost constants and transmission upgrades
capital.wind.permw<-1500000  
capital.trans.permile<-1e+06
recovery.factor <- 0.1
length.line.west_east.miles <- 200  
length.line.east_downstate.miles <- 150 
cap.line.mw <- 400  

#Calculate annual costs for wind power and each transmission line
ann.cost.wind.permw <- capital.wind.permw * recovery.factor
ann.cost.west_east.perline <- capital.trans.permile * length.line.west_east.miles * 
  recovery.factor
ann.cost.east_downstate.perline <- capital.trans.permile * length.line.east_downstate.miles * 
  recovery.factor

#Idenitify the number of variables
n.variables <- n.sites + n.regions * n.hrs + (n.regions - 1) * n.hrs + 2

#Number of constraints
n.cons <- 1 + n.regions * n.hrs + (n.regions - 1) * n.hrs + 1

#Number of bounded variables
n.bounds <- n.sites + n.regions * n.hrs + (n.regions - 1) * n.hrs + 2

#Create a for loop to evaluate different levels of wind penetration
for (j in 1:length(capacities.mw)) {
  wind.cap.total.mw <- capacities.mw[j]
  
  #Define total number of turbines
  N.t <- wind.cap.total.mw/rated.cap.per.mw
  
  #For the cost calculations, I will need the CF for Case 1
  cf.case1 <- 0.2679
  
  #Develop MILP
  
  #Create linear program model with number of rows=num of constraints and
  lp.case2 <- make.lp(nrow = n.cons, ncol = n.variables)
  
  #Specify Min
  lp.control(lp.case2, sense = "min")
  
  #Define coefficients for the objective function 
  set.objfn(lp.case2, obj = c(rep(0, n.sites), rep(1, n.regions * n.hrs), 
                              rep(0, (n.regions - 1) * n.hrs), 0, 0))
  
  #Force the number of turbines at each site and the number of new transmission lines to be an integer
  set.type(lp.case2, c((1:n.sites), (n.sites + n.regions * n.hrs + (n.regions - 
                                                                      1) * n.hrs + 1), (n.sites + n.regions * n.hrs + (n.regions - 1) * n.hrs + 
                                                                                          2)), type = "integer")
  
  #Add constraints by rows
  row.add.mode(lp.case2, "on")
  
  #Sum of number of turbines at all sites equals total number of turbines
  add.constraint(lp.case2, xt = rep(1, n.sites), type = "=", rhs = N.t, indices = c(1:n.sites))
  
  #Energy balance for West region
  for (i in 1:n.hrs) {
    add.constraint(lp.case2, xt = c(sites.wind.power.per[i, 1 + sites.west], 
                                    1, -1), type = ">=", rhs = demand.mw[i, "west"] - gen.base.west.mw, 
                   indices = c(sites.west, (n.sites + i), (n.sites + n.regions * n.hrs + 
                                                             i)))
  }
  
  #Energy balance for East region
  for (i in 1:n.hrs) {
    add.constraint(lp.case2, xt = c(sites.wind.power.per[i, 1 + sites.east], 
                                    1, 1, -1), type = ">=", rhs = demand.mw[i, "east"] - gen.base.east.mw, 
                   indices = c(sites.east, (n.sites + n.hrs + i), (n.sites + n.regions * 
                                                                     n.hrs + i), (n.sites + n.regions * n.hrs + n.hrs + i)))
  }
  
  #Energy balance for Downstate region
  for (i in 1:n.hrs) {
    add.constraint(lp.case2, xt = c(sites.wind.power.per[i, 1 + sites.downstate], 
                                    1, 1), type = ">=", rhs = demand.mw[i, "downstate"] - gen.base.downstate.mw, 
                   indices = c(sites.downstate, (n.sites + 2 * n.hrs + i), (n.sites + 
                                                                              n.regions * n.hrs + n.hrs + i)))
  }
  
  #Transmission constraints with additional lines - West-East transmission constraint
  for (i in 1:n.hrs) {
    add.constraint(lp.case2, xt = c(1, -cap.line.mw), type = "<=", rhs = t.west_east.limit.mw, 
                   indices = c((n.sites + n.regions * n.hrs + i), (n.sites + n.regions * 
                                                                     n.hrs + (n.regions - 1) * n.hrs + 1)))
  }
  
  #Transmission constraints with additional lines - East-Downstate transmission constraint
  for (i in 1:n.hrs) {
    add.constraint(lp.case2, xt = c(1, -cap.line.mw), type = "<=", rhs = t.east_downstate.limit.mw, 
                   indices = c((n.sites + n.regions * n.hrs + n.hrs + i), (n.sites + 
                                                                             n.regions * n.hrs + (n.regions - 1) * n.hrs + 2)))
  }
  
  #Cost constraint
  add.constraint(lp.case2, xt = c(rep(1, n.regions * n.hrs), (8760 * cf.case1 * 
                                                                ann.cost.west_east.perline/ann.cost.wind.permw), (8760 * cf.case1 * 
                                                                                                                    ann.cost.east_downstate.perline/ann.cost.wind.permw)), type = "<=", 
                 rhs = (sum(demand.mw$total.mw) - n.hrs * gen.base.total.mw - n.hrs * 
                          cf.case1 * wind.cap.total.mw), indices = c(((n.sites + 1):(n.sites + 
                                                                                       n.regions * n.hrs)), (n.sites + n.regions * n.hrs + (n.regions - 
                                                                                                                                              1) * n.hrs + 1), (n.sites + n.regions * n.hrs + (n.regions - 1) * 
                                                                                                                                                                  n.hrs + 2)))
  
  #Turn off row entry mode
  row.add.mode(lp.case2, "off")
  
  #Set bounds on decision variables
  set.bounds(lp.case2, lower = rep(0, n.variables), upper = c(wind.sites$n.t.max, 
                                                              rep(cap.gen.other.west.mw, n.hrs), rep(cap.gen.other.east.mw, n.hrs), 
                                                              rep(cap.gen.other.downstate.mw, n.hrs), rep(Inf, n.hrs), rep(Inf, n.hrs), 
                                                              Inf, Inf))
  
  #Solve optimization problem
  solve(lp.case2)
  
  #Get outputs
  turbine.distr <- get.variables(lp.case2)[1:n.sites]
  othergen.total.mwh <- get.objective(lp.case2)
  t.west_east.total.mwh <- sum(get.variables(lp.case2)[(n.sites + n.regions * 
                                                          n.hrs + 1):(n.sites + n.regions * n.hrs + n.hrs)])
  t.east_downstate.total.mwh <- sum(get.variables(lp.case2)[(n.sites + n.regions * 
                                                               n.hrs + n.hrs + 1):(n.sites + n.regions * n.hrs + 2 * n.hrs)])
  n.lines.west_east <- get.variables(lp.case2)[(n.sites + n.regions * n.hrs + 
                                                  (n.regions - 1) * n.hrs + 1)]
  n.lines.east_downstate <- get.variables(lp.case2)[(n.sites + n.regions * 
                                                       n.hrs + (n.regions - 1) * n.hrs + 2)]
  cf.wind <- (sum(demand.mw$total.mw) - othergen.total.mwh - gen.base.total.mw * 
                n.hrs)/(wind.cap.total.mw * 8760)
  
  #Save the turbine distribution, othergen.total and transmission totals
  assign(sprintf("turbine.distr.case2.%s", wind.cap.total.mw), turbine.distr)
  assign(sprintf("othergen.total.mwh.case2.%s", wind.cap.total.mw), othergen.total.mwh)
  assign(sprintf("t.west_east.total.mwh.case2.%s", wind.cap.total.mw), t.west_east.total.mwh)
  assign(sprintf("t.east_downstate.total.mwh.case2.%s", wind.cap.total.mw), 
         t.east_downstate.total.mwh)
  assign(sprintf("n.lines.west_east.case2.%s", wind.cap.total.mw), n.lines.west_east)
  assign(sprintf("n.lines.east_downstate.case2.%s", wind.cap.total.mw), n.lines.east_downstate)
  assign(sprintf("cf.wind.case2.%s", wind.cap.total.mw), cf.wind)
  
  #Close For loop
}

#Add capacity factors to results.cf data frame
results.cf$case2 <- cf.wind.case2.15000

#Add wind turbine distribution for 15 GW total capacity
results.turbine.distr$n.t.case2.15gw <- c(turbine.distr.case1.15000)