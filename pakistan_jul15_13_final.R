# Model of ACF in Pakistan
# David Dowdy
# Revision date: July 15, 2013

############################################################################
# DEFINE USER PARAMETERS
# Entries & TB transmission
extrapars <-c(TP = 100000,     # total population size
              B_a = 0.0001608,   # transmissibility of TB
              Bx_smn = 0.22,  # relative infectivity of smear-neg cases
              LP = 0.5)        # protection against reinfection if latently infected

# Initial population -> make sure this agrees w/ extrapars
total_pop <- 100000    # size of initial population
latent_prot <- 0.5     # protection of latent TB against reinfection
transmit <- 0.0001608
trans_smn <-0.22

# TB progression
lpfast <- 0.07         # fast progression of latent TB
fast_dur <- 2          # duration of fast progression
rel_dur <- 5 	     # duration of "recent tx" compartment
lpslow <- 0.0005       # slow progression of latent TB
sp_prop <- 0.65        # proportion of pulm TB that is smear-positive (Fig 1, exclude "not done")
ep_prop <- 0.24        # based on 2011 Indus hospital graph (slide 24) + Table (adult cases)

# Mortality
# 65-year life expectancy in Pakistan, but stop at age 50 given population structure
# start at 10 yrs old, based on pediatric data from Indus
m_tb_sp <- 0.233       # mortality from TB, smear-pos (fit to TB mort)
m_tb_sn <- m_tb_sp*0.286   # mortality from TB, smear-neg

# Self Cure
cure_sp <-0.1          # self-cure rate, smear-pos
cure_sn <-cure_sp*2.7  # self-cure, smear-neg

# TB Diagnosis & Treatment
# Pakistan prevalence/incidence ratio calculated at 364/231 = 1.576
# multiply everything by 1/0.75 to account for the 25% of population that is <10yo

# fitting to the following:
# based on adding in projected smear-positive proportion to those w/o smear at Indus,
# then adding to smear-positive vs. all-form notification rate for rest of intervention area
# assuming that extrapulmonary notification rate was constant
# 2008-2010: 54.25% smear-positive notification rate, 24% extrapulm, 21.75% smear-neg pulm
# 2011: 50.66% smear-positive, 24% extrapulm, 25.34% smear-neg pulm
# So, notification rates per 100,000:
# 2008 public: 59.3 smear-pos, 26.2 extrapulm, 23.8 smear-neg pulm (109.3 total)
# 2009 public: 74.5 smear-pos, 33.0 extrapulm, 29.9 smear-neg pulm (137.4 total)
# 2010 public: 95.7 smear-pos, 42.4 extrapulm, 38.4 smear-neg pulm (176.5 total)
# 2011 public: 173.7 smear-pos, 82.3 extrapulm, 86.9 smear-neg pulm (342.9 total)

# private sector: estimated 45% of all TB cases (Wells IJTLD 2011), 20% PPM notified (WHO 2010)
# thus, in 2008, presume that 59.3*0.8 = 55% of all cases (sm-pos)
# all cases = 59.3*0.8/0.55 = 86.3
# unreg. private sector = (86.3-59.3) = 27.0 = 0.455 times the rate measured
# at baseline, will assume that this is entirely absorbed by the intervention, from 2008-2011
#   in proportional fashion
# use 50% treatment success in private clinics as baseline: Lonnroth IJTLD 2003; 7:165
# urban may equal rural? Hoa prev survey in Vietnam
# two scenarios: equal incidence to Pakistan and 1.5x incidence

dx_pub_sp_0 <- 0.407
dx_pub_sn_0 <- 0.1585
dx_pub_ep_0 <- 0.502

prv_size <- 0.455

dx_prv_sp_0 <- prv_size*dx_pub_sp_0
dx_prv_sn_0 <- prv_size*dx_pub_sn_0
dx_prv_ep_0 <- prv_size*dx_pub_ep_0

# dx_pub_sp_m <- 0
# dx_pub_sn_m <- 0
# dx_pub_ep_m <- 0
# dx_prv_sp_m <- 0
# dx_prv_sn_m <- 0
# dx_prv_ep_m <- 0

# dx_pub_sp <- dx_pub_sp_0
# dx_pub_sn <- dx_pub_sn_0
# dx_pub_ep <- dx_pub_ep_0
# dx_prv_sp <- dx_prv_sp_0
# dx_prv_sn <- dx_prv_sn_0
# dx_prv_ep <- dx_prv_ep_0

# chg_pak <- rep(0,3)


tx_pub_0 <- 0.815217  # 0.75/0.92, removing deaths
tx_pub_m <- -0.0674   # linear regression from 2008 to 2010
                      # 2008: 0.75/0.92, 2009: 0.71/0.94, 2010: 0.66/0.97
tx_prv <- 0.5 	    # unregulated private sector tx success rate

# Relapse Rate
rel <- 0.0243	    # rate of relapse in the "recent tx" compartment


####################################################################
# Set initial matrices and population size
# Reminder of the Matrix Structure:
# matrix =7 x 9
# rows: 7 TB states (1=S, 2=Lf, 3=Ls, 4=Asp, 5=Asn, 6=Aep, 7=Recent Tx)
# cols: 9 age categories (0-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49) 

# Start all the various matrices off as empty:
statemat <- array(total_pop*0.1, dim=c(7,9)) # matrix of initial population sizes (10% per agecat)
mortmat <- array(0, dim=c(7,9))       # matrix of mortality rates
entrymat <- array(0, dim=c(7,9))      # matrix of entries into the population
mastermat <- array(0, dim=c(7,9,7,9)) # master transmission matrix, 1st 2 dimensions: loss, next 2:gain


# Set the initial population states:
#statemat[1,] <- statemat[1,]*0.75
#statemat[2,] <- statemat[2,]*(0.05-0.00364)
#statemat[3,] <- statemat[3,]*0.19
#statemat[4,] <- statemat[4,]*0.00364*sp_prop*(1-ep_prop)
#statemat[5,] <- statemat[5,]*0.00364*(1-sp_prop)*(1-ep_prop)
#statemat[6,] <- statemat[6,]*0.00364*ep_prop
#statemat[7,] <- statemat[7,]*0.01
#statemat[,1] <- statemat[,1]*2  # 1st agecat is twice as large


# now, prepare them as a list for entry into ode:
statename<-paste("S",1:63,sep="")
statevect <-c(statemat)
names(statevect)<-statename

################################################
# Mortality

mortmat[4,] <- mortmat[4,] - m_tb_sp
mortmat[5,] <- mortmat[5,] - m_tb_sn
mortmat[6,] <- mortmat[6,] - m_tb_sn
mortmat[,9] <- mortmat[,9] - 0.2

mortname<-paste("m",1:63,sep="")
mortvect <-c(mortmat)
names(mortvect)<-mortname

i <- rep(1:9,each=7)
j <- rep(1:7,9)
k <- array(c(j,i),dim=c(63,2))
m <- array(c(k,k),dim=c(63,4))
mastermat[m]<-mastermat[m]+mortvect

################################################
# Aging
n <- array(c(c(1:7),rep(1,7),c(1:7),rep(2,7)), dim=c(7,4))
mastermat[n] <- 0.1
o <- array(c(rep((1:7),7),rep(2:8,each=7),rep((1:7),7),rep(3:9,each=7)), dim=c(49,4))
mastermat[o] <- 0.2


################################################
# TB Infection & Reinfection
# TB infection occurs entirely within the diff eq function.

################################################
# TB Progression (Latent/fast-> Latent/slow-> Active)

p <- array(c(rep(2,8),c(2:9),rep(4,8),c(2:9)), dim=c(8,4))
q <- array(c(rep(2,8),c(2:9),rep(5,8),c(2:9)), dim=c(8,4))
r <- array(c(rep(2,8),c(2:9),rep(6,8),c(2:9)), dim=c(8,4))

p2 <- array(c(rep(3,8),c(2:9),rep(4,8),c(2:9)), dim=c(8,4))
q2 <- array(c(rep(3,8),c(2:9),rep(5,8),c(2:9)), dim=c(8,4))
r2 <- array(c(rep(3,8),c(2:9),rep(6,8),c(2:9)), dim=c(8,4))

s2 <- array(c(rep(2,9),c(1:9),rep(3,9),c(1:9)), dim=c(9,4))

p2b <- array(c(rep(7,8),c(2:9),rep(4,8),c(2:9)), dim=c(8,4))
q2b <- array(c(rep(7,8),c(2:9),rep(5,8),c(2:9)), dim=c(8,4))
r2b <- array(c(rep(7,8),c(2:9),rep(6,8),c(2:9)), dim=c(8,4))
s2b <- array(c(rep(7,9),c(1:9),rep(3,9),c(1:9)), dim=c(9,4))

mastermat[p]<-lpfast*sp_prop*(1-ep_prop)
mastermat[q]<-lpfast*(1-sp_prop)*(1-ep_prop)
mastermat[r]<-lpfast*ep_prop

mastermat[p2]<-lpslow*sp_prop*(1-ep_prop)
mastermat[q2]<-lpslow*(1-sp_prop)*(1-ep_prop)
mastermat[r2]<-lpslow*ep_prop

mastermat[s2]<-1/fast_dur

mastermat[p2b]<-rel*sp_prop*(1-ep_prop)
mastermat[q2b]<-rel*(1-sp_prop)*(1-ep_prop)
mastermat[r2b]<-rel*ep_prop 

mastermat[s2b]<-1/rel_dur

################################################
# Treatment

p3 <- array(c(rep(4,8),c(2:9),rep(7,8),c(2:9)), dim=c(8,4))
q3 <- array(c(rep(5,8),c(2:9),rep(7,8),c(2:9)), dim=c(8,4))
r3 <- array(c(rep(6,8),c(2:9),rep(7,8),c(2:9)), dim=c(8,4))

mastermat[p3]<-(dx_pub_sp_0*tx_pub_0) + (dx_prv_sp_0*tx_prv) + cure_sp
mastermat[q3]<-(dx_pub_sn_0*tx_pub_0) + (dx_prv_sn_0*tx_prv) + cure_sn
mastermat[r3]<-(dx_pub_ep_0*tx_pub_0) + (dx_prv_ep_0*tx_prv) + cure_sn
                                    

###################################################################
# Function for input into the differential equation solver:
# Reminder: the "c" function reads down columns of a matrix
# Thus, item [1,4] is much further down than [4,1]

# first, put mastermat into a useable format
mastername<-paste("M",1:3969,sep="")  # 63*63
mastervect <-c(mastermat)
names(mastervect)<-mastername
mastervect<-c(mastervect,extrapars)

TBdx<-function(t, state, parameters) {
        x<-as.list(c(state,parameters))
        with(x,{
             statev <- unlist(x[1:63])
             masterv <- unlist(x[64:(64+63*63-1)])

             # force of infection
             force <- B_a*((S4+S11+S18+S25+S32+S39+S46+S53+S60)+
               Bx_smn*(S5+S12+S19+S26+S33+S40+S47+S54+S61))
            
             # now, create an "infection matrix"
             infectmat <- array(0, dim=c(7,9,7,9)) 
             p4 <- array(c(rep(1,9),c(1:9),rep(2,9),c(1:9)), dim=c(9,4))
             q4 <- array(c(rep(3,9),c(1:9),rep(2,9),c(1:9)), dim=c(9,4))
		 r4 <- array(c(rep(7,9),c(1:9),rep(2,9),c(1:9)), dim=c(9,4))
             infectmat[p4]<-force
             infectmat[q4]<-force*(1-LP)
		 infectmat[r4]<-force*(1-LP)

             # convert infectmat into a vector
             infectv<-c(infectmat)

             # add infectv to masterv
             masterv3<-masterv+infectv

             # create the change matrix
             mastermatr <-array(masterv3,dim=c(63,63))

             # sum deaths, so entries can be added
             i5<-array(c(1:63),dim=c(63,2))
             i5a<-c(1:63)
             deaths<-sum(mastermatr[i5]*statev[i5a])
             # increase the population by 3%/yr
             mastermatr[i5]<-mastermatr[i5]+0.03

             # rate of change
             delta <- rep(0,63)
             delta[1] <- sum(mastermatr[,1]*statev) - 
                         sum(mastermatr[1,]*statev[1]) + 
                         mastermatr[1,1]*statev[1] -
                         deaths
             delta[2] <- sum(mastermatr[,2]*statev) - 
                         sum(mastermatr[2,]*statev[2]) + 
                         mastermatr[2,2]*statev[2]
             delta[3] <- sum(mastermatr[,3]*statev) - 
                         sum(mastermatr[3,]*statev[3]) + 
                         mastermatr[3,3]*statev[3]
             delta[4] <- sum(mastermatr[,4]*statev) - 
                         sum(mastermatr[4,]*statev[4]) + 
                         mastermatr[4,4]*statev[4]
             delta[5] <- sum(mastermatr[,5]*statev) - 
                         sum(mastermatr[5,]*statev[5]) + 
                         mastermatr[5,5]*statev[5]
             delta[6] <- sum(mastermatr[,6]*statev) - 
                         sum(mastermatr[6,]*statev[6]) + 
                         mastermatr[6,6]*statev[6] 
             delta[7] <- sum(mastermatr[,7]*statev) - 
                         sum(mastermatr[7,]*statev[7]) + 
                         mastermatr[7,7]*statev[7]
             delta[8] <- sum(mastermatr[,8]*statev) - 
                         sum(mastermatr[8,]*statev[8]) + 
                         mastermatr[8,8]*statev[8]
             delta[9] <- sum(mastermatr[,9]*statev) - 
                         sum(mastermatr[9,]*statev[9]) + 
                         mastermatr[9,9]*statev[9]
             delta[10] <- sum(mastermatr[,10]*statev) - 
                          sum(mastermatr[10,]*statev[10]) + 
                          mastermatr[10,10]*statev[10]
             delta[11] <- sum(mastermatr[,11]*statev) - 
               sum(mastermatr[11,]*statev[11]) + 
               mastermatr[11,11]*statev[11]
             delta[12] <- sum(mastermatr[,12]*statev) - 
               sum(mastermatr[12,]*statev[12]) + 
               mastermatr[12,12]*statev[12]
             delta[13] <- sum(mastermatr[,13]*statev) - 
               sum(mastermatr[13,]*statev[13]) + 
               mastermatr[13,13]*statev[13]
             delta[14] <- sum(mastermatr[,14]*statev) - 
               sum(mastermatr[14,]*statev[14]) + 
               mastermatr[14,14]*statev[14]
             delta[15] <- sum(mastermatr[,15]*statev) - 
               sum(mastermatr[15,]*statev[15]) + 
               mastermatr[15,15]*statev[15]
             delta[16] <- sum(mastermatr[,16]*statev) - 
               sum(mastermatr[16,]*statev[16]) + 
               mastermatr[16,16]*statev[16] 
             delta[17] <- sum(mastermatr[,17]*statev) - 
               sum(mastermatr[17,]*statev[17]) + 
               mastermatr[17,17]*statev[17]
             delta[18] <- sum(mastermatr[,18]*statev) - 
               sum(mastermatr[18,]*statev[18]) + 
               mastermatr[18,18]*statev[18]
             delta[19] <- sum(mastermatr[,19]*statev) - 
               sum(mastermatr[19,]*statev[19]) + 
               mastermatr[19,19]*statev[19]
             delta[20] <- sum(mastermatr[,20]*statev) - 
               sum(mastermatr[20,]*statev[20]) + 
               mastermatr[20,20]*statev[20]
             delta[21] <- sum(mastermatr[,21]*statev) - 
               sum(mastermatr[21,]*statev[21]) + 
               mastermatr[21,21]*statev[21]
             delta[22] <- sum(mastermatr[,22]*statev) - 
               sum(mastermatr[22,]*statev[22]) + 
               mastermatr[22,22]*statev[22]
             delta[23] <- sum(mastermatr[,23]*statev) - 
               sum(mastermatr[23,]*statev[23]) + 
               mastermatr[23,23]*statev[23]
             delta[24] <- sum(mastermatr[,24]*statev) - 
               sum(mastermatr[24,]*statev[24]) + 
               mastermatr[24,24]*statev[24]
             delta[25] <- sum(mastermatr[,25]*statev) - 
               sum(mastermatr[25,]*statev[25]) + 
               mastermatr[25,25]*statev[25]
             delta[26] <- sum(mastermatr[,26]*statev) - 
               sum(mastermatr[26,]*statev[26]) + 
               mastermatr[26,26]*statev[26] 
             delta[27] <- sum(mastermatr[,27]*statev) - 
               sum(mastermatr[27,]*statev[27]) + 
               mastermatr[27,27]*statev[27]
             delta[28] <- sum(mastermatr[,28]*statev) - 
               sum(mastermatr[28,]*statev[28]) + 
               mastermatr[28,28]*statev[28]
             delta[29] <- sum(mastermatr[,29]*statev) - 
               sum(mastermatr[29,]*statev[29]) + 
               mastermatr[29,29]*statev[29]
             delta[30] <- sum(mastermatr[,30]*statev) - 
               sum(mastermatr[30,]*statev[30]) + 
               mastermatr[30,30]*statev[30]
             delta[31] <- sum(mastermatr[,31]*statev) - 
               sum(mastermatr[31,]*statev[31]) + 
               mastermatr[31,31]*statev[31]
             delta[32] <- sum(mastermatr[,32]*statev) - 
               sum(mastermatr[32,]*statev[32]) + 
               mastermatr[32,32]*statev[32]
             delta[33] <- sum(mastermatr[,33]*statev) - 
               sum(mastermatr[33,]*statev[33]) + 
               mastermatr[33,33]*statev[33]
             delta[34] <- sum(mastermatr[,34]*statev) - 
               sum(mastermatr[34,]*statev[34]) + 
               mastermatr[34,34]*statev[34]
             delta[35] <- sum(mastermatr[,35]*statev) - 
               sum(mastermatr[35,]*statev[35]) + 
               mastermatr[35,35]*statev[35]
             delta[36] <- sum(mastermatr[,36]*statev) - 
               sum(mastermatr[36,]*statev[36]) + 
               mastermatr[36,36]*statev[36] 
             delta[37] <- sum(mastermatr[,37]*statev) - 
               sum(mastermatr[37,]*statev[37]) + 
               mastermatr[37,37]*statev[37]
             delta[38] <- sum(mastermatr[,38]*statev) - 
               sum(mastermatr[38,]*statev[38]) + 
               mastermatr[38,38]*statev[38]
             delta[39] <- sum(mastermatr[,39]*statev) - 
               sum(mastermatr[39,]*statev[39]) + 
               mastermatr[39,39]*statev[39]
             delta[40] <- sum(mastermatr[,40]*statev) - 
               sum(mastermatr[40,]*statev[40]) + 
               mastermatr[40,40]*statev[40]
             delta[41] <- sum(mastermatr[,41]*statev) - 
               sum(mastermatr[41,]*statev[41]) + 
               mastermatr[41,41]*statev[41]
             delta[42] <- sum(mastermatr[,42]*statev) - 
               sum(mastermatr[42,]*statev[42]) + 
               mastermatr[42,42]*statev[42]
             delta[43] <- sum(mastermatr[,43]*statev) - 
               sum(mastermatr[43,]*statev[43]) + 
               mastermatr[43,43]*statev[43]
             delta[44] <- sum(mastermatr[,44]*statev) - 
               sum(mastermatr[44,]*statev[44]) + 
               mastermatr[44,44]*statev[44]
             delta[45] <- sum(mastermatr[,45]*statev) - 
               sum(mastermatr[45,]*statev[45]) + 
               mastermatr[45,45]*statev[45]
             delta[46] <- sum(mastermatr[,46]*statev) - 
               sum(mastermatr[46,]*statev[46]) + 
               mastermatr[46,46]*statev[46] 
             delta[47] <- sum(mastermatr[,47]*statev) - 
               sum(mastermatr[47,]*statev[47]) + 
               mastermatr[47,47]*statev[47]
             delta[48] <- sum(mastermatr[,48]*statev) - 
               sum(mastermatr[48,]*statev[48]) + 
               mastermatr[48,48]*statev[48]
             delta[49] <- sum(mastermatr[,49]*statev) - 
               sum(mastermatr[49,]*statev[49]) + 
               mastermatr[49,49]*statev[49]
             delta[50] <- sum(mastermatr[,50]*statev) - 
               sum(mastermatr[50,]*statev[50]) + 
               mastermatr[50,50]*statev[50]
             delta[51] <- sum(mastermatr[,51]*statev) - 
               sum(mastermatr[51,]*statev[51]) + 
               mastermatr[51,51]*statev[51]
             delta[52] <- sum(mastermatr[,52]*statev) - 
               sum(mastermatr[52,]*statev[52]) + 
               mastermatr[52,52]*statev[52]
             delta[53] <- sum(mastermatr[,53]*statev) - 
               sum(mastermatr[53,]*statev[53]) + 
               mastermatr[53,53]*statev[53]
             delta[54] <- sum(mastermatr[,54]*statev) - 
               sum(mastermatr[54,]*statev[54]) + 
               mastermatr[54,54]*statev[54]
             delta[55] <- sum(mastermatr[,55]*statev) - 
               sum(mastermatr[55,]*statev[55]) + 
               mastermatr[55,55]*statev[55]
             delta[56] <- sum(mastermatr[,56]*statev) - 
               sum(mastermatr[56,]*statev[56]) + 
               mastermatr[56,56]*statev[56] 
             delta[57] <- sum(mastermatr[,57]*statev) - 
               sum(mastermatr[57,]*statev[57]) + 
               mastermatr[57,57]*statev[57]
             delta[58] <- sum(mastermatr[,58]*statev) - 
               sum(mastermatr[58,]*statev[58]) + 
               mastermatr[58,58]*statev[58]
             delta[59] <- sum(mastermatr[,59]*statev) - 
               sum(mastermatr[59,]*statev[59]) + 
               mastermatr[59,59]*statev[59]
             delta[60] <- sum(mastermatr[,60]*statev) - 
               sum(mastermatr[60,]*statev[60]) + 
               mastermatr[60,60]*statev[60]
             delta[61] <- sum(mastermatr[,61]*statev) - 
               sum(mastermatr[61,]*statev[61]) + 
               mastermatr[61,61]*statev[61]
             delta[62] <- sum(mastermatr[,62]*statev) - 
               sum(mastermatr[62,]*statev[62]) + 
               mastermatr[62,62]*statev[62]
             delta[63] <- sum(mastermatr[,63]*statev) - 
               sum(mastermatr[63,]*statev[63]) + 
               mastermatr[63,63]*statev[63]
             # return the rate of change
             list(c(delta))
             })
}

#############################################################
# INCIDENCE/PREVALENCE/MORTALITY CALCULATOR

outgen<-function(statevect, mastervect) {
  x<-as.list(c(statevect,mastervect))
  with(x,{
    B_a<-transmit
    statev <- unlist(x[1:63])
    masterv <- unlist(x[64:(64+63*63-1)])
    
    # force of infection
    force <- B_a*((S4+S11+S18+S25+S32+S39+S46+S53+S60)+
             Bx_smn*(S5+S12+S19+S26+S33+S40+S47+S54+S61))
    
    # now, create an "infection matrix"
    infectmat <- array(0, dim=c(7,9,7,9)) 
    p4 <- array(c(rep(1,9),c(1:9),rep(2,9),c(1:9)), dim=c(9,4))
    q4 <- array(c(rep(3,9),c(1:9),rep(2,9),c(1:9)), dim=c(9,4))
    r4 <- array(c(rep(7,9),c(1:9),rep(2,9),c(1:9)), dim=c(9,4))
    infectmat[p4]<-force
    infectmat[q4]<-force*(1-LP)
    infectmat[r4]<-force*(1-LP)
    
    # convert infectmat into a vector
    infectv<-c(infectmat)
    
    # add infectv to masterv
    masterv3<-masterv+infectv
    
    # create the change matrix
    mastermatr <-array(masterv3,dim=c(63,63))
    
    # sum deaths, so entries can be added
    i5<-array(c(1:63),dim=c(63,2))
    i5a<-c(1:63)
    deaths<-sum(mastermatr[i5]*statev[i5a])
    
    # True Incidence and Prevalence
    i2u<-c(1,8,15,22,29,36,43,50,57)
    i2lf<-i2u+1
    i2ls<-i2u+2
    i2sp<-i2u+3
    i2sn<-i2u+4
    i2ep<-i2u+5
    i2r<-i2u+6
    i2act<-c(i2sp,i2sn,i2ep)
    i2sus<-c(i2lf,i2ls,i2r)
    incarr <- array(c(rep(i2lf,3),rep(i2ls,3),rep(i2r,3),rep(i2act,3)),dim=c(81,2))
    incidence <-sum(mastermatr[incarr]*statev[incarr[,1]])
    treatarr <- array(c(i2act,rep(i2r,3)),dim=c(27,2))
    retxarr <- array(c(rep(i2r,3),i2act),dim=c(27,2))
    TBprev<-sum(statev[i2act])
    retx <- sum(mastermatr[retxarr]*statev[retxarr[,1]])

    # Notification Rate
    detect <- sum(statev[i2sp]*dx_pub_sp) + sum(statev[i2sn]*dx_pub_sn) +
              sum(statev[i2ep]*dx_pub_ep) 
    pdet_sp <- sum(statev[i2sp]*dx_pub_sp)/detect
    pdet_sn <- sum(statev[i2sn]*dx_pub_sn)/detect
    pdet_ep <- sum(statev[i2ep]*dx_pub_ep)/detect

    # TB mortality
    # TB death matrix
    deathmatr <-array(0,dim=c(63,63))
    i4a<-array(i2sp,dim=c(9,2))
    i4b<-array(i2sn,dim=c(9,2))
    i4c<-array(i2ep,dim=c(9,2))
    
    i33b<-array(i2act,dim=c(27,2))

    deathmatr[i4a]<-m_tb_sp
    deathmatr[i4b]<-m_tb_sn
    deathmatr[i4c]<-m_tb_sn

    TBmort <-sum(deathmatr[i33b]*statev[i2act])
    TBprevratio <- TBprev/incidence
    TBmortratio <- TBmort/TBprev*100
    
    i50<-c(incidence, TBprev, TBmort, retx, detect,
           pdet_sp*detect, pdet_sn*detect, pdet_ep*detect)
    i50
  })
}


##############################################################
# FIT 2009 AND 2010


times2 <- seq(0, 0.1, by = 0.1)
infile<-(read.csv("sep15_12eq.csv"))
statevect <-infile[51,]
statevect <-(as.numeric(statevect[3:65]))
statename<-paste("S",1:63,sep="")
names(statevect)<-statename
library(deSolve)
out_fin <- ode(statevect, times2, TBdx, mastervect)
out_fin2 <-array(0,dim=c(2,64))


finmat <-array(0,dim=c(71,8))
intmat <-array(0,dim=c(71,8))

dx_pub_sp_m <- 0.2205
dx_pub_sn_m <- 0.0742
dx_pub_ep_m <- 0.2922
dx_prv_sp_m <- 0
dx_prv_sn_m <- 0
dx_prv_ep_m <- 0

for (abc in 1:20) {

 statevect <-(as.numeric(out_fin[2,2:64]))
 statename<-paste("S",1:63,sep="")
 names(statevect)<-statename


dx_pub_sp <- dx_pub_sp_0 + dx_pub_sp_m * (abc/10)
dx_pub_sn <- dx_pub_sn_0 + dx_pub_sn_m * (abc/10)
dx_pub_ep <- dx_pub_ep_0 + dx_pub_ep_m * (abc/10)

dx_prv_sp <- prv_size*dx_pub_sp_0 + dx_prv_sp_m * (abc/10)
dx_prv_sn <- prv_size*dx_pub_sn_0 + dx_prv_sn_m * (abc/10)
dx_prv_ep <- prv_size*dx_pub_ep_0 + dx_prv_ep_m * (abc/10)

p3 <- array(c(rep(4,8),c(2:9),rep(7,8),c(2:9)), dim=c(8,4))
q3 <- array(c(rep(5,8),c(2:9),rep(7,8),c(2:9)), dim=c(8,4))
r3 <- array(c(rep(6,8),c(2:9),rep(7,8),c(2:9)), dim=c(8,4))

mastermat[p3]<-(dx_pub_sp*(tx_pub_0+tx_pub_m*abc/10)) + (dx_prv_sp*tx_prv) + cure_sp
mastermat[q3]<-(dx_pub_sn*(tx_pub_0+tx_pub_m*abc/10)) + (dx_prv_sn*tx_prv) + cure_sn
mastermat[r3]<-(dx_pub_ep*(tx_pub_0+tx_pub_m*abc/10)) + (dx_prv_ep*tx_prv) + cure_sn

mastername<-paste("M",1:3969,sep="")  # 63*63
mastervect <-c(mastermat)
names(mastervect)<-mastername
mastervect<-c(mastervect,extrapars)

out_fin <- ode(statevect, times2, TBdx, mastervect)
finmat[abc,]<-outgen(out_fin[2,2:64],mastervect)
intmat[abc,]<-outgen(out_fin[2,2:64],mastervect)
}

out_2010 <-out_fin[2,2:64]
out_fin2 <-out_fin

####################################
# BASELINE

int_sp <- 0
int_sn <- 0
int_ep <- 0

int_prv <- 1.0

out_fin2[2,2:64] <-out_2010

dx_pub_sp <- dx_pub_sp_0 + dx_pub_sp_m * (20/10)
dx_pub_sn <- dx_pub_sn_0 + dx_pub_sn_m * (20/10) 
dx_pub_ep <- dx_pub_ep_0 + dx_pub_ep_m * (20/10)
dx_prv_sp <- int_prv*(prv_size*dx_pub_sp_0 + dx_prv_sp_m * (20/10))
dx_prv_sn <- int_prv*(prv_size*dx_pub_sn_0 + dx_prv_sn_m * (20/10))
dx_prv_ep <- int_prv*(prv_size*dx_pub_ep_0 + dx_prv_ep_m * (20/10))

for (abc in 21:70) {

 statevect <-(as.numeric(out_fin2[2,2:64]))
 statename<-paste("S",1:63,sep="")
 names(statevect)<-statename

dx_pub_sp <- dx_pub_sp + max(0, dx_pub_sp_m * (3-abc/10)*0.1)
dx_pub_sn <- dx_pub_sn + max(0, dx_pub_sp_m * (3-abc/10)*0.1)
dx_pub_ep <- dx_pub_ep + max(0, dx_pub_sp_m * (3-abc/10)*0.1)

p3 <- array(c(rep(4,8),c(2:9),rep(7,8),c(2:9)), dim=c(8,4))
q3 <- array(c(rep(5,8),c(2:9),rep(7,8),c(2:9)), dim=c(8,4))
r3 <- array(c(rep(6,8),c(2:9),rep(7,8),c(2:9)), dim=c(8,4))

mastermat[p3]<-(dx_pub_sp*(tx_pub_0+tx_pub_m*20/10)) + (dx_prv_sp*tx_prv) + cure_sp
mastermat[q3]<-(dx_pub_sn*(tx_pub_0+tx_pub_m*20/10)) + (dx_prv_sn*tx_prv) + cure_sn
mastermat[r3]<-(dx_pub_ep*(tx_pub_0+tx_pub_m*20/10)) + (dx_prv_ep*tx_prv) + cure_sn

mastername<-paste("M",1:3969,sep="")  # 63*63
mastervect <-c(mastermat)
names(mastervect)<-mastername
mastervect<-c(mastervect,extrapars)

out_fin2 <- ode(statevect, times2, TBdx, mastervect)
finmat[abc,]<-outgen(out_fin2[2,2:64],mastervect)

}

#############################################
# INTERVENTION


int_sp <- 1.939
int_sn <- 1.2627
int_ep <- 2.042

int_prv <- 0

out_fin2[2,2:64] <-out_2010
for (abc in 21:70) {

 statevect <-(as.numeric(out_fin2[2,2:64]))
 statename<-paste("S",1:63,sep="")
 names(statevect)<-statename


dx_pub_sp <- dx_pub_sp_0 + dx_pub_sp_m * (20/10) + min(int_sp*(abc-20)/10,int_sp)
dx_pub_sn <- dx_pub_sn_0 + dx_pub_sn_m * (20/10) + min(int_sn*(abc-20)/10,int_sn)
dx_pub_ep <- dx_pub_ep_0 + dx_pub_ep_m * (20/10) + min(int_ep*(abc-20)/10,int_ep)

dx_prv_sp <- max(0,(3-abc/10)*(prv_size*dx_pub_sp_0 + dx_prv_sp_m * (20/10)))
dx_prv_sn <- max(0,(3-abc/10)*(prv_size*dx_pub_sn_0 + dx_prv_sn_m * (20/10)))
dx_prv_ep <- max(0,(3-abc/10)*(prv_size*dx_pub_ep_0 + dx_prv_ep_m * (20/10)))

p3 <- array(c(rep(4,8),c(2:9),rep(7,8),c(2:9)), dim=c(8,4))
q3 <- array(c(rep(5,8),c(2:9),rep(7,8),c(2:9)), dim=c(8,4))
r3 <- array(c(rep(6,8),c(2:9),rep(7,8),c(2:9)), dim=c(8,4))

mastermat[p3]<-(dx_pub_sp*(tx_pub_0+tx_pub_m*20/10)) + (dx_prv_sp*tx_prv) + cure_sp
mastermat[q3]<-(dx_pub_sn*(tx_pub_0+tx_pub_m*20/10)) + (dx_prv_sn*tx_prv) + cure_sn
mastermat[r3]<-(dx_pub_ep*(tx_pub_0+tx_pub_m*20/10)) + (dx_prv_ep*tx_prv) + cure_sn

mastername<-paste("M",1:3969,sep="")  # 63*63
mastervect <-c(mastermat)
names(mastervect)<-mastername
mastervect<-c(mastervect,extrapars)

out_fin2 <- ode(statevect, times2, TBdx, mastervect)
intmat[abc,]<-outgen(out_fin2[2,2:64],mastervect)

}

fit11 <-rep(0,8)
for (abcd in 1:8) {
fit11[abcd] <-sum(intmat[21:30,abcd])/10
}

c("detected (374.68)", "smpos (189.81)", "smneg (94.94)", "EP (89.92)")
fit11[5:8]

# write.csv(intmat,file="sep15_12int.csv")
# write.csv(finmat,file="sep15_12fin.csv")



