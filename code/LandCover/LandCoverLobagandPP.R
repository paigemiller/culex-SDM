r<-raster('C:/Users/Robert/OneDrive/Documents/Robbie/8910 ENM/Project/Africa_IVC_20130316_final_MG.tif')

global<-raster('C:/Users/Robert/OneDrive/Documents/Robbie/8910 ENM/Project/Global_Ecophysiography.tif')

coverraster<-raster('C:/Users/Robert/OneDrive/Documents/Robbie/8910 ENM/Project/Globcover2009_V2.3_Global_/GLOBCOVER_L4_200901_200912_V2.3.tif')

cells<-cellFromRowCol(data,-50:-55,0:2)




colnames(pipiens.coords)<-c("x",'y')


EFvalues<-extract(global, pipiens.coords)

covervalues<-extract(coverraster, pipiens.coords)
pipiensenvt<-extract(environmental.data.rs, pipiens.coords)

coverbackground<-extract(coverraster, backgr.pts)

coverbg<-matrix(0, nrow=nrow(backgr.pts), ncol=23)

cover<-matrix(0, nrow=nrow(pipiens.coords), ncol=23)

for (i in 1:length(covervalues))
{
  if (covervalues[i]==11) cover[i,1]<-1
  if (covervalues[i]==14) cover[i,2]<-1
  if (covervalues[i]==20) cover[i,3]<-1
  if (covervalues[i]==30) cover[i,4]<-1
  if (covervalues[i]==40) cover[i,5]<-1
  if (covervalues[i]==50) cover[i,6]<-1
  if (covervalues[i]==60) cover[i,7]<-1
  if (covervalues[i]==70) cover[i,8]<-1
  if (covervalues[i]==90) cover[i,9]<-1
  if (covervalues[i]==100) cover[i,10]<-1
  if (covervalues[i]==110) cover[i,11]<-1
  if (covervalues[i]==120) cover[i,12]<-1
  if (covervalues[i]==130) cover[i,13]<-1
  if (covervalues[i]==140) cover[i,14]<-1
  if (covervalues[i]==150) cover[i,15]<-1
  if (covervalues[i]==160) cover[i,16]<-1
  if (covervalues[i]==170) cover[i,17]<-1
  if (covervalues[i]==180) cover[i,18]<-1
  if (covervalues[i]==190) cover[i,19]<-1
  if (covervalues[i]==200) cover[i,20]<-1
  if (covervalues[i]==210) cover[i,21]<-1
  if (covervalues[i]==220) cover[i,22]<-1
  if (covervalues[i]==230) cover[i,23]<-1
  
}

for (i in 1:length(coverbackground))
{
  if (coverbackground[i]==11) coverbg[i,1]<-1
  if (coverbackground[i]==14) coverbg[i,2]<-1
  if (coverbackground[i]==20) coverbg[i,3]<-1
  if (coverbackground[i]==30) coverbg[i,4]<-1
  if (coverbackground[i]==40) coverbg[i,5]<-1
  if (coverbackground[i]==50) coverbg[i,6]<-1
  if (coverbackground[i]==60) coverbg[i,7]<-1
  if (coverbackground[i]==70) coverbg[i,8]<-1
  if (coverbackground[i]==90) coverbg[i,9]<-1
  if (coverbackground[i]==100) coverbg[i,10]<-1
  if (coverbackground[i]==110) coverbg[i,11]<-1
  if (coverbackground[i]==120) coverbg[i,12]<-1
  if (coverbackground[i]==130) coverbg[i,13]<-1
  if (coverbackground[i]==140) coverbg[i,14]<-1
  if (coverbackground[i]==150) coverbg[i,15]<-1
  if (coverbackground[i]==160) coverbg[i,16]<-1
  if (coverbackground[i]==170) coverbg[i,17]<-1
  if (coverbackground[i]==180) coverbg[i,18]<-1
  if (coverbackground[i]==190) coverbg[i,19]<-1
  if (coverbackground[i]==200) coverbg[i,20]<-1
  if (coverbackground[i]==210) coverbg[i,21]<-1
  if (coverbackground[i]==220) coverbg[i,22]<-1
  if (coverbackground[i]==230) coverbg[i,23]<-1
  
}

y<-c()
for(i in 1:23)
{
  if (sum(coverbg[,i])==0)
    y<-rbind(y,i)
}
coverbg<-coverbg[,-y]

background.cover<-cbind(data.background,coverbg)



x<-c()
for(i in 1:23)
{
   if (sum(cover[,i])==0)
     x<-rbind(x,i)
}
cover<-cover[,-x]

pipiensall<-cbind(pipiensenvt, cover)
pipiens.char<-data.frame(cbind(pipiensenvt, covervalues))
pipiens.char$covervalues<-as.character(pipiens.char$covervalues)


z <- which(is.na(pipiensall[,1])) 
pipiens.presenceall <- pipiensall[-z,]
n.presence <- dim(pipiens.presenceall)[1] 
train.presence.cover <- pipiens.presenceall[c(train.presence.id),] 
train.background.cover <- background.cover[c(train.background.id),] 
test.presence.cover <- pipiens.presenceall[-c(train.presence.id),] 
test.background.cover <- background.cover[-c(train.presence.id),]

#test.backgrond.cover[,87:]



z <- which(is.na(pipiens.char[,1])) 
pipiens.presence.char <- pipiens.char[-z,]
n.presence <- dim(pipiens.presence.char)[1] 
train.presence.char <- pipiens.presence.char[c(train.presence.id),] 
#train.background.cover <- data.background[c(train.background.id),] 
test.presence.char<- pipiens.presence.char[-c(train.presence.id),] 
#test.background.cover <- data.background[-c(train.presence.id),]

KDE<-pp.kde(p=train.presence, bgrd=train.background)

predictions<-predict.pp.kde(model=KDE, x=test.presence)

predictkde<-predict.pp.kde(KDE, x=environmental.data.rs)
### RANGEBAG
Rmodel<-rb(train.presence)
rastersum<- predict(environmental.data.rs,Rmodel[[1]])
for (i in 2:100)
{
  rasterpredict<-predict(environmental.data.rs,model[[i]])
  
  rastersum<-rastersum+rasterpredict
}
RBrastermean<-rastersum/length(model)
plot(RBrastermean)


###

maxmodel<-maxent(environmental.data.rs,p=pipiens.coords, a=backgr.pts)
predictionmax<-predict(environmental.data.rs,maxmodel)

plot(predictionmax)

###
model<-lobag.oc(p=train.presence.char, n.votes=100)

predictions.char<-predictSvm(model=model, p=test.presence.char)

###
model<-lobag.oc(p=train.presence.cover, n.votes=100)



predictions.cover<-predictSvm(model=model, p=test.presence.cover,a=test.background.cover)


##
model<-lobag.oc(p=train.presence, n.votes=100)
predictions<-predictSvm(model=model, p=test.presence)


rastersum<- predict(environmental.data.rs,model[[1]])
for (i in 2:100)
{
  rasterpredict<-predict(environmental.data.rs,model[[i]])
  
  rastersum<-rastersum+rasterpredict
}
rastermean<-rastersum/length(model)
plot(rastermean)


mean(predictions.cover$p.out)

mean(predictions$p.out)

#####Quinquefasciatus
covervalues<-extract(coverraster, quinq.coords)
quinqenvt<-extract(environmental.data.rs, quinq.coords)

cover.quinq<-matrix(0, nrow=nrow(quinq.coords), ncol=23)

for (i in 1:length(covervalues))
{
  if (covervalues[i]==11) cover.quinq[i,1]<-1
  if (covervalues[i]==14) cover.quinq[i,2]<-1
  if (covervalues[i]==20) cover.quinq[i,3]<-1
  if (covervalues[i]==30) cover.quinq[i,4]<-1
  if (covervalues[i]==40) cover.quinq[i,5]<-1
  if (covervalues[i]==50) cover.quinq[i,6]<-1
  if (covervalues[i]==60) cover.quinq[i,7]<-1
  if (covervalues[i]==70) cover.quinq[i,8]<-1
  if (covervalues[i]==90) cover.quinq[i,9]<-1
  if (covervalues[i]==100) cover.quinq[i,10]<-1
  if (covervalues[i]==110) cover.quinq[i,11]<-1
  if (covervalues[i]==120) cover.quinq[i,12]<-1
  if (covervalues[i]==130) cover.quinq[i,13]<-1
  if (covervalues[i]==140) cover.quinq[i,14]<-1
  if (covervalues[i]==150) cover.quinq[i,15]<-1
  if (covervalues[i]==160) cover.quinq[i,16]<-1
  if (covervalues[i]==170) cover.quinq[i,17]<-1
  if (covervalues[i]==180) cover.quinq[i,18]<-1
  if (covervalues[i]==190) cover.quinq[i,19]<-1
  if (covervalues[i]==200) cover.quinq[i,20]<-1
  if (covervalues[i]==210) cover.quinq[i,21]<-1
  if (covervalues[i]==220) cover.quinq[i,22]<-1
  if (covervalues[i]==230) cover.quinq[i,23]<-1
  
}
x<-c()
for(i in 1:23)
{
  if (sum(cover.quinq[,i])==0)
    x<-rbind(x,i)
}
cover.quinq<-cover.quinq[,-x]

quinqall<-cbind(quinqenvt, cover.quinq)

z <- which(is.na(quinqall[,1])) 
quinq.presenceall <- quinqall[-z,]
n.presence <- dim(quinq.presenceall)[1]

train.presence.id.quinq <- sample(seq(1,n.presence),ceiling(0.8*n.presence)) 
train.presence.quinq <- quinqall[c(train.presence.id.quinq),]
test.presence.quinq <- quinqall[-c(train.presence.id.quinq),]

 
train.presence.cover.quinq <- quinq.presenceall[c(train.presence.id.quinq),] 
#train.background.cover <- data.background[c(train.background.id),] 
test.presence.cover.quinq <- quinq.presenceall[-c(train.presence.id.quinq),] 
#test.background.cover <- data.background[-c(train.presence.id),]

model<-lobag.oc(p=train.presence.cover.quinq, n.votes=100)

predictions.cover.quinq<-predictSvm(model=model, p=test.presence.cover.quinq)

model<-lobag.oc(p=train.presence.cover.quinq[,1:86], n.votes=100)

predictions.quinq<-predictSvm(model=model, p=test.presence.cover.quinq[,1:86])

mean(predictions.cover.quinq$p.out)

mean(predictions.quinq$p.out)


Na.remove<-function(x){
  if (is.na(x)==T){
    x<-0
  }else{
    x<-x
  }
  
}
T1992noNA<-apply(T1992, c(1,2), FUN=Na.remove)
