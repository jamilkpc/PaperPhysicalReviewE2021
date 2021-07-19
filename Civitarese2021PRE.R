### Figure 1

library(deSolve)
ode_model = function(t,u,p){
  #p[1] = m
  #p[2] = s
  #p[3] = q
  if (u>1/2){
    dc = (1-u)*((1-p[1])*(u^p[3]) + p[1]*p[2]/2 + p[1]*(1-p[2])) - u*((1-p[1])*((1-u)^p[3]) + p[1]*p[2]/2)
  } else {
    dc = (1-u)*((1-p[1])*(u^p[3]) + p[1]*p[2]/2) - u*((1-p[1])*((1-u)^p[3]) + p[1]*p[2]/2  + p[1]*(1-p[2]))
  }
  return(list(dc))
}

optSc <- function(m,s,q,u0){
  p = c(m,s,q)
  sol = ode(u0,c(0:10000)/10,ode_model,p)
  sc = sol[10000,2]
  return(sc)
}

#par(mfrow=c(3,1), xpd=TRUE)
n = 50
qvec = c(5,9,9)
uvec = c(4,4,1)/10
tvec = c('a','b','c')
for (s in 1:3){
  matSc = matrix(0,(n-1),(n-1))
  q = qvec[s]
  u0 = uvec[s]
  for (i in 2:n){
    for (j in 2:n){
      c1 = (i-1)/n
      c2 = (j-1)/n
      matSc[(i-1),(j-1)] = optSc(c1,c2,q,u0)
    }
  }
  filled.contour(matSc, xlab = "m", ylab="s", cex.lab=1.35, cex.axis = 1.15)
  mtext(bquote('('*.(tvec[s])*') q = '*.(qvec[s])*'; '*c[0]*' = '*.(uvec[s])), line = 1, adj = 0.35, cex = 1.3)
}


### Figure 2
par(mar=c(5, 4, 0.5, 4), xpd=TRUE)
ind = c()
minM = c()
for (i in 1:20){
  ind = append(ind,i)
  q = i
  c = (0:10000)/10000
  x = min(c^q + (1-c)^q - q*c*(1-c)^(q-1) - q*(1-c)*c^(q-1))
  y = x/(x-1)
  minM = append(minM,y)
}
q = ind
plot(q,minM, type = "b", xlim = c(0,20), ylim = c(0,0.35), ylab = 'm*', cex.lab = 1.35, cex.axis = 1.05)

### Figure 3
par(mar=c(5, 4, 0.5, 8), xpd=TRUE)

library(dplyr)

m = 1/5
q = 13
cL = (0:4999)/10000
sL = 2*cL + 2*((1-m)/m)*(cL*(1-cL)^q - (1-cL)*cL^q)
cH = (5001:10000)/10000
sH = 2*(1-cH) - 2*((1-m)/m)*(cH*(1-cH)^q - (1-cH)*cH^q)
s = append(append(sL,NA),sH)
s = na_if((s<0)*10000+(s >= 0)*s,10000)
s = na_if((s>1)*10000+(s <= 1)*s,10000)
plot(cbind(s,(0:10000)/10000), type = 'l',ylim=c(0,1), xlim=c(0,1), xlab = "s", ylab = "Concentration",  col='red', lwd = 1, lty = 2, cex.lab=1.35)

m = 1/20
q = 13
cL = (0:4999)/10000
sL = 2*cL + 2*((1-m)/m)*(cL*(1-cL)^q - (1-cL)*cL^q)
cH = (5001:10000)/10000
sH = 2*(1-cH) - 2*((1-m)/m)*(cH*(1-cH)^q - (1-cH)*cH^q)
s = append(append(sL,NA),sH)
s = na_if((s<0)*10000+(s >= 0)*s,10000)
s = na_if((s>1)*10000+(s <= 1)*s,10000)
plot(cbind(s,(0:10000)/10000), type = 'l',ylim=c(0,1), xlim=c(0,1), xlab = "s", ylab = "Concentration",  col='red', lwd = 1, lty = 2, cex.lab=1.35)

q = 2
cL = (0:4999)/10000
sL = 2*cL + 2*((1-m)/m)*(cL*(1-cL)^q - (1-cL)*cL^q)
cH = (5001:10000)/10000
sH = 2*(1-cH) - 2*((1-m)/m)*(cH*(1-cH)^q - (1-cH)*cH^q)
s = append(append(sL,NA),sH)
s = na_if((s<0)*10000+(s >= 0)*s,10000)
s = na_if((s>1)*10000+(s <= 1)*s,10000)
lines(cbind(s,(0:10000)/10000), col = 'red', lwd = 1, lty = 1)

m = 1/3
q = 13
cL = (0:4999)/10000
sL = 2*cL + 2*((1-m)/m)*(cL*(1-cL)^q - (1-cL)*cL^q)
cH = (5001:10000)/10000
sH = 2*(1-cH) - 2*((1-m)/m)*(cH*(1-cH)^q - (1-cH)*cH^q)
s = append(append(sL,NA),sH)
s = na_if((s<0)*10000+(s >= 0)*s,10000)
s = na_if((s>1)*10000+(s <= 1)*s,10000)
#plot(cbind(s,(0:10000)/10000), type = 'l',ylim=c(0,1), ylab = "Concentration", xlab = "External Field Strength", col="red")
lines(cbind(s,(0:10000)/10000), col = 1, lwd = 1, lty = 2)

q = 2
cL = (0:4999)/10000
sL = 2*cL + 2*((1-m)/m)*(cL*(1-cL)^q - (1-cL)*cL^q)
cH = (5001:10000)/10000
sH = 2*(1-cH) - 2*((1-m)/m)*(cH*(1-cH)^q - (1-cH)*cH^q)
s = append(append(sL,NA),sH)
s = na_if((s<0)*10000+(s >= 0)*s,10000)
s = na_if((s>1)*10000+(s <= 1)*s,10000)
lines(cbind(s,(0:10000)/10000), col = 1, lwd = 1, lty = 1)

legend('topright',inset=c(-0.25,0.3), c("q = 2", "q = 13"), lty = c(2,1), col = c('red','red'), bty = 'n',
       title = "m = 1/20", cex = 1.1)
legend('topright',inset=c(-0.25,0.6), c("q = 2", "q = 13"), lty = c(2,1), col = c('black','black'), bty = 'n',
       title = "m = 1/3", cex = 1.1)

### Figure 4
par(mar=c(5, 4.5, 0.5, 8), mfrow=c(2,1), xpd=TRUE)

s = 0.5
q = 2
cL = (0:4999)/10000
mL = (-(1-cL)*cL^q + cL*(1-cL)^q)/(-(1-cL)*cL^q + cL*(1-cL)^q - cL + s/2)
cH = (5001:10000)/10000
mH = (-(1-cH)*cH^q + cH*(1-cH)^q)/(-(1-cH)*cH^q + cH*(1-cH)^q + 1 - cH - s/2)
m = append(append(mL,NA),mH)
m = na_if((m<0)*10000+(m >= 0)*m,10000)
m = na_if((m>1)*10000+(m <= 1)*m,10000)
plot(cbind(m,(0:10000)/10000), type = 'l',ylim=c(0,1), ylab = "Concentration", xlab = "m", lty = 1, col="orange", cex.lab = 1.6, cex.axis = 1.1)

s = 0.8
q = 2
cL = (0:4999)/10000
mL = (-(1-cL)*cL^q + cL*(1-cL)^q)/(-(1-cL)*cL^q + cL*(1-cL)^q - cL + s/2)
cH = (5001:10000)/10000
mH = (-(1-cH)*cH^q + cH*(1-cH)^q)/(-(1-cH)*cH^q + cH*(1-cH)^q + 1 - cH - s/2)
m = append(append(mL,NA),mH)
m = na_if((m<0)*10000+(m >= 0)*m,10000)
m = na_if((m>1)*10000+(m <= 1)*m,10000)
lines(cbind(m,(0:10000)/10000), lty = 1, col = "red")

s = 0.995
q = 2
cL = (0:4999)/10000
mL = (-(1-cL)*cL^q + cL*(1-cL)^q)/(-(1-cL)*cL^q + cL*(1-cL)^q - cL + s/2)
cH = (5001:10000)/10000
mH = (-(1-cH)*cH^q + cH*(1-cH)^q)/(-(1-cH)*cH^q + cH*(1-cH)^q + 1 - cH - s/2)
m = append(append(mL,NA),mH)
m = na_if((m<0)*10000+(m >= 0)*m,10000)
m = na_if((m>1)*10000+(m <= 1)*m,10000)
lines(cbind(m,(0:10000)/10000), lty = 1)

legend('right',inset=c(-0.225,0), c("s = 0.5", "s = 0.8", "s = 0.995"), lty = c(1,1,1), col = c('orange','red','black'), bty = 'n', title = 'q = 2', cex = 1.2)

s = 0.5
q = 13
cL = (0:4999)/10000
mL = (-(1-cL)*cL^q + cL*(1-cL)^q)/(-(1-cL)*cL^q + cL*(1-cL)^q - cL + s/2)
cH = (5001:10000)/10000
mH = (-(1-cH)*cH^q + cH*(1-cH)^q)/(-(1-cH)*cH^q + cH*(1-cH)^q + 1 - cH - s/2)
m = append(append(mL,NA),mH)
m = na_if((m<0)*10000+(m >= 0)*m,10000)
m = na_if((m>1)*10000+(m <= 1)*m,10000)
plot(cbind(m,(0:10000)/10000), type = 'l',ylim=c(0,1), ylab = "Concentration", xlab = "m", lty = 1, col="orange", cex.lab = 1.6, cex.axis = 1.1)

s = 0.8
q = 13
cL = (0:4999)/10000
mL = (-(1-cL)*cL^q + cL*(1-cL)^q)/(-(1-cL)*cL^q + cL*(1-cL)^q - cL + s/2)
cH = (5001:10000)/10000
mH = (-(1-cH)*cH^q + cH*(1-cH)^q)/(-(1-cH)*cH^q + cH*(1-cH)^q + 1 - cH - s/2)
m = append(append(mL,NA),mH)
m = na_if((m<0)*10000+(m >= 0)*m,10000)
m = na_if((m>1.00001)*10000+(m <= 1.00001)*m,10000)
lines(cbind(m,(0:10000)/10000), lty = 1 ,col = "red")

s = 0.995
q = 13
cL = (0:4999)/10000
mL = (-(1-cL)*cL^q + cL*(1-cL)^q)/(-(1-cL)*cL^q + cL*(1-cL)^q - cL + s/2)
cH = (5001:10000)/10000
mH = (-(1-cH)*cH^q + cH*(1-cH)^q)/(-(1-cH)*cH^q + cH*(1-cH)^q + 1 - cH - s/2)
m = append(append(mL,NA),mH)
m = na_if((m<0)*10000+(m >= 0)*m,10000)
m = na_if((m>1.00001)*10000+(m <= 1.00001)*m,10000)
lines(cbind(m,(0:10000)/10000), lty = 1)

legend('right',inset=c(-0.225,0), c("s = 0.5", "s = 0.8", "s = 0.995"), lty = c(1,1,1), col = c('orange','red','black'), bty = 'n', title = 'q = 13', cex = 1.2)
dev.off()


### Figure 5
par(mfrow=c(2,2))

s = 0.9
q = 5
cL = (0:4999)/10000
mL = (-(1-cL)*cL^q + cL*(1-cL)^q)/(-(1-cL)*cL^q + cL*(1-cL)^q - cL + s/2)
cH = (5001:10000)/10000
mH = (-(1-cH)*cH^q + cH*(1-cH)^q)/(-(1-cH)*cH^q + cH*(1-cH)^q + 1 - cH - s/2)
m = append(append(mL,NA),mH)
m = na_if((m<0)*10000+(m >= 0)*m,10000)
m = na_if((m>1)*10000+(m <= 1)*m,10000)
plot(cbind(m,(0:10000)/10000), col = "black", type = "l", ylab = "Concentration", xlab = "External Field Strength", cex.lab = 1.35, cex.axis = 1.1)
mtext('(a) q = 5 and s = 0.9', line = 1, cex = 1.2)
lines((0:10000)/10000,rep(0.5,10001), col = "red", lty = 2)
arrows(0.0, 0.49, 0, 0.01, length = 0.1, col = "grey")
arrows(0.0, 0.51, 0, 0.99, length = 0.1, col = "grey")

arrows(0.175, 0.49, 0.175, 0.15, length = 0.1, col = "grey")
arrows(0.175, 0.51, 0.175, 0.85, length = 0.1, col = "grey")
arrows(0.175, 0.00, 0.175, 0.125, length = 0.1, col = "grey")
arrows(0.175, 1.00, 0.175, 0.875, length = 0.1, col = "grey")

arrows(0.35, 0.49, 0.35, 0.42, length = 0.1, col = "grey")
arrows(0.35, 0.51, 0.35, 0.58, length = 0.1, col = "grey")
arrows(0.35, 0.00, 0.35, 0.40, length = 0.1, col = "grey")
arrows(0.35, 1.00, 0.35, 0.60, length = 0.1, col = "grey")

arrows(0.55, 0.49, 0.55, 0.45, col = "grey", length = 0.1)
arrows(0.55, 0.51, 0.55, 0.55, col = "grey", length = 0.1)
arrows(0.55, 0, 0.55, 0.43, col = "grey", length = 0.1)
arrows(0.55, 1, 0.55, 0.57, col = "grey", length = 0.1)

arrows(0.75, 0.49, 0.75, 0.46, col = "grey", length = 0.1)
arrows(0.75, 0.51, 0.75, 0.54, col = "grey", length = 0.1)
arrows(0.75, 0, 0.75, 0.44, col = "grey", length = 0.1)
arrows(0.75, 1, 0.75, 0.56, col = "grey", length = 0.1)

arrows(0.95, 0.49, 0.95, 0.46, col = "grey", length = 0.1)
arrows(0.95, 0.51, 0.95, 0.54, col = "grey", length = 0.1)
arrows(0.95, 0, 0.95, 0.44, col = "grey", length = 0.1)
arrows(0.95, 1, 0.95, 0.56, col = "grey", length = 0.1)

s = 0.9
q = 9
cL = (0:1383)/10000
mL = (-(1-cL)*cL^q + cL*(1-cL)^q)/(-(1-cL)*cL^q + cL*(1-cL)^q - cL + s/2)
c1 = (1384:3660)/1000
m1 = rep(NA,length(c1))
cN = (3661:4999)/10000
mN = (-(1-cN)*cN^q + cN*(1-cN)^q)/(-(1-cN)*cN^q + cN*(1-cN)^q - cN + s/2)
mD = append(append(mL,m1),mN)

cN = (5001:6339)/10000
mN = (-(1-cN)*cN^q + cN*(1-cN)^q)/(-(1-cN)*cN^q + cN*(1-cN)^q + 1 - cN - s/2)
m1 = rep(NA,length(c1))
cH = (8617:10000)/10000
mH = (-(1-cH)*cH^q + cH*(1-cH)^q)/(-(1-cH)*cH^q + cH*(1-cH)^q + 1 - cH - s/2)
mU = append(append(mN,m1),mH)

m = append(append(mD,NA),mU)
m = na_if((m<0)*10000+(m >= 0)*m,10000)
m = na_if((m>1)*10000+(m <= 1)*m,10000)

plot(cbind(m,(0:10000)/10000), col = "black", type = "l", ylab = "Concentration", xlab = "External Field Strength", cex.lab = 1.35, cex.axis = 1.1)
mtext('(b) q = 9 and s = 0.9', line = 1, cex = 1.2)

cL = (0:1383)/10000
mL = rep(NA,length(cL))
c1 = (1384:3660)/10000
m1 = (-(1-c1)*c1^q + c1*(1-c1)^q)/(-(1-c1)*c1^q + c1*(1-c1)^q - c1 + s/2)
cN = (3661:4999)/10000
mN = rep(NA,length(cN))
mD = append(append(mL,m1),mN)

cN = (5001:6339)/10000
mN = rep(NA,length(cN))
c1 = (6340:8616)/10000
m1 = (-(1-c1)*c1^q + c1*(1-c1)^q)/(-(1-c1)*c1^q + c1*(1-c1)^q + 1 - c1 - s/2)
cH = (8617:10000)/10000
mH = rep(NA,length(cH))
mU = append(append(mN,m1),mH)

m = append(append(mD,NA),mU)
m = na_if((m<0)*10000+(m >= 0)*m,10000)
m = na_if((m>1)*10000+(m <= 1)*m,10000)

lines(cbind(m,(0:10000)/10000), col = "black", type = "l", lty = 2)
lines((0:10000)/10000,rep(0.5,10001), col = "red", lty = 2)
arrows(0.0, 0.49, 0.0, 0.01, col = "grey", length= 0.1)
arrows(0.0, 0.51, 0.0, 0.99, col = "grey", length= 0.1)

arrows(0.075, 0.49, 0.075, 0.42, col = "grey", length = 0.1)
arrows(0.075, 0.51, 0.075, 0.58, col = "grey", length = 0.1)
arrows(0.075, 0.27, 0.075, 0.06, col = "grey", length = 0.1)
arrows(0.075, 0.73, 0.075, 0.94, col = "grey", length = 0.1)
arrows(0.075, 0.33, 0.075, 0.395, col = "grey", length = 0.1)
arrows(0.075, 0.67, 0.075, 0.605, col = "grey", length = 0.1)
arrows(0.075, 0.00, 0.075, 0.03, col = "grey", length = 0.1)
arrows(0.075, 1.00, 0.075, 0.97, col = "grey", length = 0.1)

arrows(0.175, 0.49, 0.175, 0.45, col = "grey", length = 0.1)
arrows(0.175, 0.51, 0.175, 0.55, col = "grey", length = 0.1)
arrows(0.175, 0, 0.175, 0.43, col = "grey", length = 0.1)
arrows(0.175, 1, 0.175, 0.57, col = "grey", length = 0.1)

arrows(0.35, 0.49, 0.35, 0.46, col = "grey", length = 0.1)
arrows(0.35, 0.51, 0.35, 0.54, col = "grey", length = 0.1)
arrows(0.35, 0, 0.35, 0.44, col = "grey", length = 0.1)
arrows(0.35, 1, 0.35, 0.56, col = "grey", length = 0.1)

arrows(0.55, 0.49, 0.55, 0.46, col = "grey", length = 0.1)
arrows(0.55, 0.51, 0.55, 0.54, col = "grey", length = 0.1)
arrows(0.55, 0, 0.55, 0.44, col = "grey", length = 0.1)
arrows(0.55, 1, 0.55, 0.56, col = "grey", length = 0.1)

arrows(0.75, 0.49, 0.75, 0.46, col = "grey", length = 0.1)
arrows(0.75, 0.51, 0.75, 0.54, col = "grey", length = 0.1)
arrows(0.75, 0, 0.75, 0.44, col = "grey", length = 0.1)
arrows(0.75, 1, 0.75, 0.56, col = "grey", length = 0.1)

arrows(0.95, 0.49, 0.95, 0.46, col = "grey", length = 0.1)
arrows(0.95, 0.51, 0.95, 0.54, col = "grey", length = 0.1)
arrows(0.95, 0, 0.95, 0.44, col = "grey", length = 0.1)
arrows(0.95, 1, 0.95, 0.56, col = "grey", length = 0.1)

s = 0.5
q = 5
cL = (0:4999)/10000
mL = (-(1-cL)*cL^q + cL*(1-cL)^q)/(-(1-cL)*cL^q + cL*(1-cL)^q - cL + s/2)
cH = (5001:10000)/10000
mH = (-(1-cH)*cH^q + cH*(1-cH)^q)/(-(1-cH)*cH^q + cH*(1-cH)^q + 1 - cH - s/2)
m = append(append(mL,NA),mH)
m = na_if((m<0)*10000+(m >= 0)*m,10000)
m = na_if((m>1)*10000+(m <= 1)*m,10000)
plot(cbind(m,(0:10000)/10000), col = "black", type = "l", ylab = "Concentration", xlab = "External Field Strength", cex.lab = 1.35, cex.axis = 1.1)
mtext('(c) q = 5 and s = 0.5', line = 1, cex = 1.2)
lines((0:10000)/10000,rep(0.5,10001), col = "red", lty = 2)

arrows(0.0, 0.49, 0.0, 0.01, col = "grey", length= 0.1)
arrows(0.0, 0.51, 0.0, 0.99, col = "grey", length= 0.1)

arrows(0.175, 0.49, 0.175, 0.07, col = "grey", length= 0.1)
arrows(0.175, 0.51, 0.175, 0.93, col = "grey", length= 0.1)
arrows(0.175, 0.00, 0.175, 0.05, col = "grey", length= 0.1)
arrows(0.175, 1.00, 0.175, 0.95, col = "grey", length= 0.1)

arrows(0.35, 0.49, 0.35, 0.145, col = "grey", length= 0.1)
arrows(0.35, 0.51, 0.35, 0.855, col = "grey", length= 0.1)
arrows(0.35, 0.00, 0.35, 0.125, col = "grey", length= 0.1)
arrows(0.35, 1.00, 0.35, 0.875, col = "grey", length= 0.1)

arrows(0.55, 0.49, 0.55, 0.21, col = "grey", length= 0.1)
arrows(0.55, 0.51, 0.55, 0.79, col = "grey", length= 0.1)
arrows(0.55, 0.00, 0.55, 0.19, col = "grey", length= 0.1)
arrows(0.55, 1.00, 0.55, 0.81, col = "grey", length= 0.1)

arrows(0.75, 0.49, 0.75, 0.24, col = "grey", length= 0.1)
arrows(0.75, 0.51, 0.75, 0.76, col = "grey", length= 0.1)
arrows(0.75, 0.00, 0.75, 0.22, col = "grey", length= 0.1)
arrows(0.75, 1.00, 0.75, 0.78, col = "grey", length= 0.1)

arrows(0.95, 0.49, 0.95, 0.26, col = "grey", length= 0.1)
arrows(0.95, 0.51, 0.95, 0.74, col = "grey", length= 0.1)
arrows(0.95, 0.00, 0.95, 0.24, col = "grey", length= 0.1)
arrows(0.95, 1.00, 0.95, 0.76, col = "grey", length= 0.1)

s = 0.5
q = 9
cL = (0:4999)/10000
mL = (-(1-cL)*cL^q + cL*(1-cL)^q)/(-(1-cL)*cL^q + cL*(1-cL)^q - cL + s/2)
cH = (5001:10000)/10000
mH = (-(1-cH)*cH^q + cH*(1-cH)^q)/(-(1-cH)*cH^q + cH*(1-cH)^q + 1 - cH - s/2)
m = append(append(mL,NA),mH)
m = na_if((m<0)*10000+(m >= 0)*m,10000)
m = na_if((m>1)*10000+(m <= 1)*m,10000)
plot(cbind(m,(0:10000)/10000), col = "black", type = "l", ylab = "Concentration", xlab = "External Field Strength", cex.lab = 1.35, cex.axis = 1.1)
mtext('(d) q = 9 and s = 0.5', line = 1, cex = 1.2)
lines((0:10000)/10000,rep(0.5,10001), col = "red", lty = 2)

arrows(0.0, 0.49, 0.0, 0.01, col = "grey", length= 0.1)
arrows(0.0, 0.51, 0.0, 0.99, col = "grey", length= 0.1)

arrows(0.175, 0.49, 0.175, 0.09, col = "grey", length= 0.1)
arrows(0.175, 0.51, 0.175, 0.91, col = "grey", length= 0.1)
arrows(0.175, 0.00, 0.175, 0.07, col = "grey", length= 0.1)
arrows(0.175, 1.00, 0.175, 0.93, col = "grey", length= 0.1)

arrows(0.35, 0.49, 0.35, 0.21, col = "grey", length= 0.1)
arrows(0.35, 0.51, 0.35, 0.79, col = "grey", length= 0.1)
arrows(0.35, 0.00, 0.35, 0.19, col = "grey", length= 0.1)
arrows(0.35, 1.00, 0.35, 0.81, col = "grey", length= 0.1)

arrows(0.55, 0.49, 0.55, 0.24, col = "grey", length= 0.1)
arrows(0.55, 0.51, 0.55, 0.76, col = "grey", length= 0.1)
arrows(0.55, 0.00, 0.55, 0.22, col = "grey", length= 0.1)
arrows(0.55, 1.00, 0.55, 0.78, col = "grey", length= 0.1)

arrows(0.75, 0.49, 0.75, 0.25, col = "grey", length= 0.1)
arrows(0.75, 0.51, 0.75, 0.75, col = "grey", length= 0.1)
arrows(0.75, 0.00, 0.75, 0.23, col = "grey", length= 0.1)
arrows(0.75, 1.00, 0.75, 0.77, col = "grey", length= 0.1)

arrows(0.95, 0.49, 0.95, 0.26, col = "grey", length= 0.1)
arrows(0.95, 0.51, 0.95, 0.74, col = "grey", length= 0.1)
arrows(0.95, 0.00, 0.95, 0.24, col = "grey", length= 0.1)
arrows(0.95, 1.00, 0.95, 0.76, col = "grey", length= 0.1)


# Figure 6
o_traj = function(v_init){
  media_strength = v_init[1]
  proportion_independent = v_init[2]
  initial_division = v_init[3]
  q = v_init[4]
  steps = v_init[5]
  size = 200
  independence_par = 0.5
  
  opinion = rbinom(size,1,initial_division)
  sc = rep(NA,(size+1))
  
  sc[1] = sum(opinion)/size
  p = rbinom(steps,1,independence_par)
  s = rbinom(steps,1,proportion_independent)
  for (i in 1:steps){
    noise = sample(size,(q+1))
    if (rbinom(1,1,media_strength)==0){
      opinion[noise[1]] = if(sum(opinion[noise[2:(q+1)]])%%q==0){opinion[noise[2]]} else{opinion[noise[1]]}
      sc[(i+1)] = sum(opinion)/size
    } else {
      opinion[noise[1]] = (1-s[i])*floor(median(opinion)) + s[i]*p[i]
      opinion[noise[1]] = (1-s[i])*opinion[noise[1]] + s[i]*p[i]
      sc[(i+1)] = sum(opinion)/size
    }
  }
  return(sc)
}

par(mar = c(4, 5, 4, 4), mfrow = c(2,1), xpd=TRUE)
#a1 = o_traj(c(0.08,0.8,0.10,11,5000000))
plot((1:5000001)/200, a1, type = "l", col = "red", ylab = "Concentration", xlab = "Monte Carlo Steps", ylim=c(0,1), cex.axis = 1.4, cex.lab = 1.85)
grid()
mtext(bquote('(a) m = 0.08; q = 11; s = 0.80'), line = 1, cex = 2)

#a2 = o_traj(c(0.16,0.8,0.10,11,5000000))
plot((1:5000001)/200, a2, type = "l", col = "red", ylab = "Concentration", xlab = "Monte Carlo Steps", ylim=c(0,1), cex.axis = 1.4, cex.lab = 1.85)
grid()
mtext(bquote('(b) m = 0.16; q = 11; s = 0.80'), line = 1, cex = 2)

#### Figure 7
o_count = function(v_init){
  media_strength = v_init[1]
  proportion_independent = v_init[2]
  initial_division = v_init[3]
  q = v_init[4]
  size = 200
  threshold = 200/3
  count = 0
  independence_par = 0.5
  
  opinion = rbinom(size,1,initial_division)
  sto = sum(opinion)
  
  steps = 5000
  p = rbinom(steps,1,independence_par)
  s = rbinom(steps,1,proportion_independent)
  
  for (i in 1:steps){
    noise = sample(size,(q+1))
    if (rbinom(1,1,media_strength)==0){
      opinion[noise[1]] = if(sum(opinion[noise[2:(q+1)]])%%q==0){opinion[noise[2]]} else{opinion[noise[1]]}
      if ((100-abs(100-sum(opinion)))>threshold){count=i}
    } else {
      opinion[noise[1]] = (1-s[i])*floor(median(opinion)) + s[i]*p[i]
      if ((100-abs(100-sum(opinion)))>threshold){count=i}
    }
  }
  return(count)
}

o_avgcount = function(v_init){
  opinion = apply(t(matrix(rep(v_init,200),4,200)),1,o_count)
  return(median(opinion))
}

library(future.apply)

paramvec = cbind(((0:10)/20),rep(0.1,11),rep(0.5,11),rep(5,11))
sto = future_apply(paramvec, 1, o_avgcount, future.seed = TRUE)
sto

paramvec = cbind(((0:10)/20),rep(0.1,11),rep(0.5,11),rep(9,11))
sto2 = future_apply(paramvec, 1, o_avgcount, future.seed = TRUE)
sto2

paramvec = cbind(((0:10)/20),rep(0.3,11),rep(0.5,11),rep(5,11))
sto3 = future_apply(paramvec, 1, o_avgcount, future.seed = TRUE)
sto3

paramvec = cbind(((0:10)/20),rep(0.3,11),rep(0.5,11),rep(9,11))
sto4 = future_apply(paramvec, 1 , o_avgcount, future.seed = TRUE)
sto4

plot((0:10)/20,log(sto), type = "l", ylim=c(4,9), ylab = "Log of Time to Supermajority", xlab = "m", cex.lab = 1.25)
lines((0:10)/20,log(sto2), type = "l", lty = 2)
lines((0:10)/20,log(sto3), type = "l", col = "blue")
lines((0:10)/20,log(sto4), type = "l", col = "blue", lty = 2)
lines((0:10)/20,rep(log(5000),11), col = "chartreuse4", lty = 2)