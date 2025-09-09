## @knitr all.functions
seg.func<-function(inter.set)
{#inter.set<-inter.DS
  fn<-length(inter.set);
  fseg<-NULL;k.beg<-k.end<-inter.set[1]; ct<-1
  for(i in 2:(fn))
    if ((inter.set[i]-inter.set[i-1])>1) {
      k.end<-inter.set[i-1];
      while((k.end-k.beg)>11)
      {fseg[[ct]]<-k.beg:(k.beg+11);
      ct<-ct+1;
      k.beg<-(k.beg+12)
      }
      fseg[[ct]]<-k.beg:k.end;
      ct<-ct+1;
      k.beg<-inter.set[i];
    }
  k.end<-inter.set[fn];
  while((k.end-k.beg)>11)
  {fseg[[ct]]<-k.beg:(k.beg+11);
  ct<-ct+1;
  k.beg<-(k.beg+12)
  }
  fseg[[ct]]<-k.beg:k.end;
  len.temp<-sapply(fseg, length); 
  if(any(len.temp<=2)) fseg<-fseg[-which(len.temp<=2)]
  return(fseg)
}

comp.func<-function(fsam) #fsam is a list of subsamples ordered with time.
  # this function is to calculate the freq of <  
  # fsam<-dat.temp.x;
{
  fn<-length(fsam)
  return(lapply(2:fn, function(x) outer(fsam[[x-1]], fsam[[x]], '<')))
}
### 
hpxy2.func<-function(outer.mat)
{ # outer.mat<-outer.x
  temp.mat<-lapply(outer.mat, function(x) t(x)%*%x/nrow(x)/ncol(x)/(ncol(x)-1))
  return(sapply(temp.mat, function(x) sum(x)-sum(diag(x))))
}
###
hp2xy.func<-function(outer.mat)
{ # outer.mat<-outer.x
  temp.mat<-lapply(outer.mat, function(x) x%*%t(x)/nrow(x)/ncol(x)/(nrow(x)-1))
  return(sapply(temp.mat, function(x) sum(x)-sum(diag(x))))
}
###
hpxy12.func<-function(outer.mat) 
{
  fn<-length(outer.mat)
  return(sapply(2:fn, function(x) 
    sum(outer.mat[[x-1]]%*%outer.mat[[x]])/nrow(outer.mat[[x-1]])/ncol(outer.mat[[x]])/nrow(outer.mat[[x]])))
}

test.stat<-function(list.x,list.y)
{
  outer.x<-comp.func(list.x); outer.y<-comp.func(list.y);
  hpx<-sapply(outer.x, mean); hpy<-sapply(outer.y, mean); 
  hpx2<-hpxy2.func(outer.x); hpy2<-hpxy2.func(outer.y); 
  hp2x<-hp2xy.func(outer.x); hp2y<-hp2xy.func(outer.y); 
  hsx1<-hpx2-hpx^2; hsx2<-hp2x-hpx^2;
  hsy1<-hpy2-hpy^2; hsy2<-hp2y-hpy^2;
  hpx12<-hpxy12.func(outer.x);   hsx12<-hpx12-hpx[-1]*hpx[-length(hpx)];      
  hpy12<-hpxy12.func(outer.y);   hsy12<-hpy12-hpy[-1]*hpy[-length(hpy)]
  
  nx<-sapply(list.x, length); ny<-sapply(list.y, length); 
  hcx<-nx/sum(nx); hcy<-ny/sum(ny); 
  len.seg<-length(hcx);
  
  hlam<-sum(nx)/(sum(nx)+sum(ny));
  cell.index.up<-cbind(1:(len.seg-2), 2:(len.seg-1));
  cell.index.down<-cbind( 2:(len.seg-1), 1:(len.seg-2));
  
  hSig.x<-diag(hsx1/hcx[-len.seg]+hsx2/hcx[-1])
  hSig.x[cell.index.up]<- hSig.x[cell.index.down] <-hsx12/hcx[-c(1,len.seg)]
  hSig.y<-diag(hsy1/hcy[-len.seg]+hsy2/hcy[-1])
  hSig.y[cell.index.up]<- hSig.y[cell.index.down] <-hsy12/hcy[-c(1,len.seg)]
  
  Ox<-sapply(outer.x, sum); Oy<-sapply(outer.y, sum)
  nl2<-sapply(outer.x, function(x) nrow(x)*ncol(x));
  ml2<-sapply(outer.y, function(y) nrow(y)*ncol(y));
  Nl<-nl2+ml2; Rxl<-nl2/Nl;
  
  hpl<-(Ox+Oy)/Nl;
  hD<-diag(hlam*sqrt(hcx[-len.seg]*hcx[-1]*(1-Rxl)/hpl/(1-hpl)))
  if(any(hD==Inf)) hD[which(hD==Inf)]<-0
  
  Exl<-(Ox+Oy)*Rxl;Eyl<-(Ox+Oy)*(1-Rxl);
  M.temp<-(Ox-Exl)^2*(1/Exl+1/(nl2-Exl))+(Oy-Eyl)^2*(1/Eyl+1/(ml2-Eyl));
  if(any(is.nan(M.temp))) M.temp[is.nan(M.temp)]<-0;
  M.stat<-sum(M.temp);
  test.stat<-M.stat/(sum(nx+ny));
  
  hsig<-hD%*%(hSig.x/hlam+hSig.y/(1-hlam))%*%hD
  eig.lam<-eigen(hsig)$values
  if(any(eig.lam<=0)) {eig.lam<-eig.lam[-which(eig.lam<=0)]; print('negative eigent value, check!')}
  
  p.val<-CompQuadForm::farebrother(test.stat, eig.lam)$Qq
  return(c(test.stat,p.val))
}

####### Analysis starts here

#fn<-file.choose();

## @knitr raw.data
fn<-"/Users/wangy/Documents/Research /Research 2025/Real Estate project/Real_Estate_Sales_2001-2021_GL.csv"
raw.dat<-read.csv(fn);

library("tidyverse")
dat1<-raw.dat[,1:10];
sum(complete.cases(dat1)); dim(dat1); 
towns<-c('Bridgeport', 'Bristol', 'Danbury', 'Hartford', 'New Haven', 
         'Norwalk', 'Stamford',  'Waterbury', 'West Hartford')
dat2<-dat1%>%filter(Town%in%towns);
dat2$Date<-as.Date(dat2$Date.Recorded, "%m/%d/%Y")
saveRDS(dat1,file='dat1.rds')
saveRDS(dat2,file='dat2.rds')
############################################# 
## @knitr read.dat1
dat1<-readRDS(file='./processed.dat/dat1.rds')
summary(dat1);

## @knitr read.dat2
dat2<-readRDS(file='./processed.dat/dat2.rds') ## this dataset contains records of bigger towns in CT
summary(dat2);

## @knitr date.process
#dat2$Date<-as.Date(dat2$Date.Recorded, "%m/%d/%Y")
min.time<-min(dat2$Date); max.time<-max(dat2$Date);
#as.numeric(difftime(dat2$Date[1:10], min.time, units = "month"))
mon.mat<-(strsplit(dat2$Date.Recorded, '/')%>%unlist%>%as.numeric%>%matrix(ncol=3, byrow=TRUE))[,-2]
dat2$mon.2001.9<-(mon.mat[,2]-2001)*12+(mon.mat[,1]-9);

hist(dat2$mon.count)
summary(dat2$mon.count);

table(dat2$Residential.Type,dat2$Property.Type)
dat3<-dat2[,c(3,4,7,9,12)];
saveRDS(dat3,file='dat3.rds')

dat3<-readRDS(file='./processed.dat/dat3.rds');

towname.i<-'Bristol'
##
dat.temp<-dat2[dat2$Town==towname.i,c(7,9)]
## @knitr prepare.dat34
aggregate(Sale.Amount~Property.Type,data=dat3, FUN = function(x) c(count=length(x),Mean = mean(x, na.rm = TRUE), Median = median(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE)))

flag.4<-!dat3$Property.Type%in%c('Commercial','Industrial', 'Public Utility','Vacant Land')
dat4<-dat3[flag.4,]
aggregate(Sale.Amount~Property.Type,data=dat4, FUN = function(x) c(count=length(x),Mean = mean(x, na.rm = TRUE), Median = median(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE)))
saveRDS(dat4,file='dat4.rds')

## @knitr dat4
dat4<-readRDS(file='./processed.dat/dat4.rds');

tab.property.time.town<-aggregate(Property.Type~mon.2001.9+Town,data=dat4, FUN = function(x) count=length(x))
colnames(tab.property.time.town)<-c('Months','Town','Sales')
ggplot2::ggplot(tab.property.time.town,                            
       aes(x = Months,
           y = Sales,
           col = Town)) + geom_line()

sum(tab.property.time.town[,3])

boxplot(Sales~Town, data=tab.property.time.town)

## @knitr toward analysis and simulation
dat4<-readRDS(file='dat4.rds');
# Let us pick Danbury and Stamford, since both are muc above zeros compared with others
dat4.sel<-dat4[,c('Town','Sale.Amount','mon.2001.9')]


dat5.s1<-split(dat4.sel, dat4.sel$Town);
dat5.s2<-lapply(dat5.s1, FUN=function(x) split(x[,2:3],x$mon.2001.9) )

dat4.DS<-dat4[dat4$Town %in% c('Danbury','Stamford'),];

par(mfrow=c(2,1))
boxplot(Sale.Amount~mon.2001.9, data=dat4.DS[dat4.DS$Town=='Danbury',], 
        ylim=c(0,10^6), main='Danbury Sales', xlab='Months since Sep 2001')
boxplot(Sale.Amount~mon.2001.9, data=dat4.DS[dat4.DS$Town=='Stamford',], 
        ylim=c(0,10^6), main='Stamford Sales', xlab='Months since Sep 2001')
par(mfrow=c(1,1))

p <- ggplot(dat4.DS, aes(x=as.factor(mon.2001.9), y=Sale.Amount,color=Town)) 
p + geom_boxplot(position = position_identity())+coord_cartesian(ylim = c(0,2.5*10^6))


dat5.s3<-lapply(dat5.s2, FUN=function(x) lapply(x, FUN=function(y) y$Sale.Amount))

saveRDS(dat5.s3,file='dat5.rds')

## @knitr dat5
dat5<-readRDS(file='./processed.dat/dat5.rds');
#towns<-c('Bridgeport', 'Bristol', 'Danbury', 'Hartford', 'New Haven', 'Norwalk',
#         'Stamford',  'Waterbury', 'West Hartford')
#two.city<-c('Danbury','Stamford');
#two.city<-c('Bridgeport','Waterbury');

two.city<-c('Bridgeport','Hartford');
## @knitr dat5.DS.gen
dat5.DS<-dat5[two.city]
dat4.DS<-dat4[dat4$Town %in% two.city,];

inter.DS<-intersect(sort(unique(dat4.DS[dat4.DS$Town==two.city[1],]$mon.2001.9)),
                    sort(unique(dat4.DS[dat4.DS$Town==two.city[2],]$mon.2001.9)));

seg.inter.DS<-seg.func(inter.DS)

seq.id<-1;

## @knitr loop.bgn
for(seq.id in c(1:length(seg.inter.DS))) { print(seg.inter.DS[[seq.id]])

## @knitr dat5.DS.v1
seq.temp<-seg.inter.DS[[seq.id]];
# seq.temp<-205:(205+21)
dat4.DS.temp<-dat4.DS[dat4.DS$mon.2001.9 %in% seq.temp,]

dat5.DS.temp.s1<-split(dat4.DS.temp, dat4.DS.temp$Town);
dat5.DS.temp.s2<-lapply(dat5.DS.temp.s1, FUN=function(x) split(x[,2:3],x$mon.2001.9) )

## @knitr dat5.DS.v1.plot
par(mfrow=c(2,1))
boxplot(Sale.Amount~mon.2001.9, data=dat4.DS.temp[dat4.DS.temp$Town==two.city[1],], 
        ylim=c(0,10^6), main=paste(two.city[1],'Sales'), xlab='Months since Sep 2001')
boxplot(Sale.Amount~mon.2001.9, data=dat4.DS.temp[dat4.DS.temp$Town==two.city[2],], 
        ylim=c(0,10^6), main=paste(two.city[2],'Sales'), xlab='Months since Sep 2001')
par(mfrow=c(1,1))

## @knitr dat5.DS.v2
p <- ggplot(dat4.DS.temp, aes(x=as.factor(mon.2001.9), y=Sale.Amount,color=Town)) 
p + geom_boxplot()+coord_cartesian(ylim = c(0,2.5*10^6))+xlab("Months since Sep 2001")

## @knitr test.stat.func
dat5.DS.temp.s3<-lapply(dat5.DS.temp.s2, FUN=function(x) lapply(x, FUN=function(y) y$Sale.Amount))
dat.temp.x<-dat5.DS.temp.s3[[1]];dat.temp.y<-dat5.DS.temp.s3[[2]];

res<-round(test.stat(dat.temp.x, dat.temp.y),3)
#print(res)
## @knitr loop.end
}


## @knitr calculation result
#(11,(12) for two.city<-c('Danbury','Stamford');

#1,2,8,14 for two.city<-c('Hartford','New Haven');
# 3/10/11 two.city<-c('Bridgeport','Stamford');






## @knitr details

outer.x<-comp.func(dat.temp.x); outer.y<-comp.func(dat.temp.y);

hpx<-sapply(outer.x, mean); hpy<-sapply(outer.y, mean); 


hpx2<-hpxy2.func(outer.x); hpy2<-hpxy2.func(outer.y); 
hp2x<-hp2xy.func(outer.x); hp2y<-hp2xy.func(outer.y); 

hsx1<-hpx2-hpx^2; hsx2<-hp2x-hpx^2;
hsy1<-hpy2-hpy^2; hsy2<-hp2y-hpy^2;


hpx12<-hpxy12.func(outer.x);   hsx12<-hpx12-hpx[-1]*hpx[-length(hpx)];      
hpy12<-hpxy12.func(outer.y);   hsy12<-hpy12-hpy[-1]*hpy[-length(hpy)]
 

nx<-sapply(dat.temp.x, length); ny<-sapply(dat.temp.y, length); 
hcx<-nx/sum(nx); hcy<-ny/sum(ny); 
len.seg<-length(hcx);

hlam<-sum(nx)/(sum(nx)+sum(ny));
cell.index.up<-cbind(1:(len.seg-2), 2:(len.seg-1));
cell.index.down<-cbind( 2:(len.seg-1), 1:(len.seg-2));

hSig.x<-diag(hsx1/hcx[-len.seg]+hsx2/hcx[-1])
hSig.x[cell.index.up]<- hSig.x[cell.index.down] <-hsx12/hcx[-c(1,len.seg)]


hSig.y<-diag(hsy1/hcy[-len.seg]+hsy2/hcy[-1])
hSig.y[cell.index.up]<- hSig.y[cell.index.down] <-hsy12/hcy[-c(1,len.seg)]


## @knitr test.stats

Ox<-sapply(outer.x, sum); Oy<-sapply(outer.y, sum)
nl2<-sapply(outer.x, function(x) nrow(x)*ncol(x));
ml2<-sapply(outer.y, function(y) nrow(y)*ncol(y));
Nl<-nl2+ml2; Rxl<-nl2/Nl;

hpl<-(Ox+Oy)/Nl;
hD<-diag(hlam*sqrt(hcx[-len.seg]*hcx[-1]*(1-Rxl)/hpl/(1-hpl)))

Exl<-(Ox+Oy)*Rxl;Eyl<-(Ox+Oy)*(1-Rxl)
M.stat<-sum((Ox-Exl)^2*(1/Exl+1/(nl2-Exl))+(Oy-Eyl)^2*(1/Eyl+1/(ml2-Eyl)));
test.stat<-M.stat/(sum(nx+ny));


hsig<-hD%*%(hSig.x/hlam+hSig.y/(1-hlam))%*%hD
eig.lam<-eigen(hsig)$values

p.val<-farebrother(test.stat, eig.lam)$Qq
print(c(test.stat,p.val))
















