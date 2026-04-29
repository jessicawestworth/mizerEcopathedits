#Finding boundary points
p_grouped<-readRDS("/Users/jessicawestworth/Desktop/BH simulations/p_grouped.rds")

#remove all c = 0, as this is just changing effort of normal fishing, and all alpha 0 as this is just status quo fishing
o<-p_grouped%>%filter(c!=0, alpha !=0)

#Y
Y_greater<-o%>%filter(BH_Y>status_Y)
#greater for every BH simulation

#B
B_greater<-o%>%filter(BH_B>status_B)
x<-c(0.2,0.3,0.4,0.5,0.6,0.9)
y<-c(1,0.6,0.4,0.3,0.2,0.1)

m1<-lm(log(y)~log(x))
summary(m1)

plot(y~x)
y1<-exp(-2.37509) * (x^-1.52954)
points(y1~x, type="l")

#SSB
SSB_greater<-o%>%filter(BH_SSB>status_SSB)

#record every point where alpha is maxed for given c value, with no repeating alphas
x<-c(0.1,0.2,0.3,0.4)
y<-c(1.0,0.4,0.2,0.1)

m1<-lm(log(y)~log(x))
summary(m1)

plot(y~x)
y2<-exp(-3.654) * (x^-1.623)
points(y2~x, type="l")

#N
N_greater<-o%>%filter(BH_N>status_N)

x<-c(0.6,0.7,0.8,0.9,2,3,4)
y<-c(1,0.9,0.8,0.7,0.3,0.2,0.1)

m1<-lm(log(y)~log(x))
summary(m1)

plot(y~x)
y3<-exp(-0.50267) * (x^-1.15159)
points(y3~x, type="l")

#LFY
LFY_greater<-o%>%filter(BH_LFY>status_LFY)
#greater for no BH simulations

#LFB
LFB_greater<-o%>%filter(BH_LFB>status_LFB)
#greater for all BH simulations

#plot boundaries on same graph
x<-seq(from=10, to =0, length.out =100)
y1<-exp(-2.37509) * (x^-1.52954)
y2<-exp(-3.654) * (x^-1.623)
y3<-exp(-0.50267) * (x^-1.15159)
plot(y1~x, type="n", ylim=c(0.1,1), xlab="c", ylab="alpha")
points(y1~x, type="l")
points(y2~x, type="l", col="blue")
points(y3~x, type="l", col="red")
#maximum yield total
o%>%filter(BH_Y==max(o$BH_Y))
#maximum yield BH fully weighted
full<-o%>%filter(alpha==1)
full%>%filter(BH_Y==max(full$BH_Y))

