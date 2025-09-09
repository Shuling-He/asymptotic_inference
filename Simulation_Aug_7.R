################# functions ##################
## @knitr all.functions
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


############ simulated data ######
#nl<-100; ml<-2*nl; k<-4;
#gam.val<-c(9,7,4,6);


#list.x<-lapply(c(1:k), function(x) rexp(nl, gam.val[x]));
#list.y<-lapply(c(1:k), function(y) rexp(ml, gam.val[y]));

simu.exp<-function(N, nl=nl)
{
  sapply(1:N, 
         function(v) {list.x<-lapply(c(1:k), function(x) rexp(nl, gam.val.x[x]));
                      list.y<-lapply(c(1:k), function(y) rexp(ml, gam.val.y[y]));
                      temp<-test.stat(list.x,list.y);
                      if ((v%%1000)==0) cat(v);
                      return(temp)
                      })
}


simu.norm<-function(N, nl=nl)
{
  sapply(1:N, 
         function(v) {list.x<-lapply(c(1:k), function(x) rnorm(nl, mu.val.x[x],sd.val.x[x]));
                      list.y<-lapply(c(1:k), function(y) rnorm(ml, mu.val.y[y],sd.val.y[y]));
                       temp<-test.stat(list.x,list.y);
                       if ((v%%1000)==0) cat(v);
                       return(temp)
         })
}

## @knitr all.simu

### simulation part 5 exp_4_diff ################
for(nl in c(20,50,100,200))
{
  print(paste ('Sample size', nl))
  rep.N<-10^4;set.seed(20250808);
  ml<-2*nl; gam.val.x<-c(9,7,4,6); k<-length(gam.val.x);
  gam.val.y<-c(9,8,4,6);
  out1<-simu.exp(rep.N, nl=nl);
  if(any(out1[2,]>1)) print('Negative Eignevalue')
  print(mean(out1[2,]<=0.05))
  file.name<-paste('Exp_diff_4_', nl,'.rds',sep='')
  saveRDS(out1, file=file.name)
}
### simulation part 6 exp_8_diff ################
rep.N<-10^4;
for(nl in c(20,50,100,200))
{
  print(paste ('Sample size', nl))
  set.seed(20250808);
  ml<-2*nl;  gam.val.x<-c(9,7,4,6,3,2,5,8); k<-length(gam.val.x);
  gam.val.y<-c(9,8,4,6,3,2,5,8);
  out1<-simu.exp(rep.N, nl=nl);
  if(any(out1[2,]>1)) print('Negative Eignevalue')
  print(mean(out1[2,]<=0.05))
  file.name<-paste('Exp_diff_8_', nl,'.rds',sep='')
  saveRDS(out1, file=file.name)
}
### simulation part 7 norm_low_diff_mean ################
for(nl in c(20,50,100,200))
{
  print(paste ('Sample size', nl))
  rep.N<-10^4;set.seed(20250808);
  mu.val.x<-c(9,7,4,6);  ml<-2*nl;  k<-length(mu.val.x);
  mu.val.y<-c(9,8,4,6);
  sd.val.x<-sd.val.y<-rep(5,k)
  out1<-simu.norm(rep.N, nl=nl);
  if(any(out1[2,]>1)) print('Negative Eignevalue')
  print(mean(out1[2,]<=0.05))
  file.name<-paste('Norm_4_lowvar_', nl,'.rds',sep='')
  saveRDS(out1, file=file.name)
}

### simulation part 8 norm_high_diff_mean ################
for(nl in c(20,50,100,200))
{
  print(paste ('Sample size', nl))
  rep.N<-10^4;set.seed(20250808);
  mu.val.x<-c(9,7,4,6); ml<-2*nl;  k<-length(mu.val.x);
  mu.val.y<-c(9,8,4,6);
  sd.val.x<-sd.val.y<-rep(10,k)
  out1<-simu.norm(rep.N, nl=nl);
  if(any(out1[2,]>1)) print('Negative Eignevalue')
  print(mean(out1[2,]<=0.05))
  file.name<-paste('Norm_4_highvar_', nl,'.rds',sep='')
  saveRDS(out1, file=file.name)
}
### simulation part 9 norm_low_diff_var ################
for(nl in c(20,50,100,200))
{
  print(paste ('Sample size', nl))
  rep.N<-10^4;set.seed(20250808);
  mu.val.x<-mu.val.y<-c(9,7,4,6);  ml<-2*nl;  k<-length(mu.val.x);
  sd.val.x<-rep(5,k)
  sd.val.y<-c(5,10,5,5)
  out1<-simu.norm(rep.N, nl=nl);
  if(any(out1[2,]>1)) print('Negative Eignevalue')
  print(mean(out1[2,]<=0.05))
  file.name<-paste('Norm_4_diffvar_', nl,'.rds',sep='')
  saveRDS(out1, file=file.name)
}

### simulation part 1 exp_4 ################
for(nl in c(20,50,100,200,500))
{
  print(paste ('Sample size', nl))
  rep.N<-10^4;set.seed(20250808);
  ml<-2*nl; gam.val.x<-gam.val.y<-c(9,7,4,6); k<-length(gam.val.x);
  out1<-simu.exp(rep.N, nl=nl);
  if(any(out1[2,]>1)) print('Negative Eignevalue')
  print(mean(out1[2,]<=0.05))
  file.name<-paste('Exp_same_', nl,'.rds',sep='')
  saveRDS(out1, file=file.name)
}

### simulation part2 for exp dist 8 ################
rep.N<-10^4;
for(nl in c(20,50,100,200))
{
  print(paste ('Sample size', nl))
  set.seed(20250808);
  ml<-2*nl;  gam.val.x<-gam.val.y<-c(9,7,4,6,3,2,5,8); k<-length(gam.val.x);
  out1<-simu.exp(rep.N, nl=nl);
  if(any(out1[2,]>1)) print('Negative Eignevalue')
  print(mean(out1[2,]<=0.05))
  file.name<-paste('Exp_same_8_', nl,'.rds',sep='')
  saveRDS(out1, file=file.name)
}


### simulation part 3 Normal short ################
for(nl in c(20,50,100,200))
{
  print(paste ('Sample size', nl))
  rep.N<-10^4;set.seed(20250808);
  mu.val.x<-mu.val.y<-c(9,7,4,6);  ml<-2*nl;  k<-length(mu.val.x);
  sd.val.x<-sd.val.y<-rep(5,k)
  out1<-simu.norm(rep.N, nl=nl);
  if(any(out1[2,]>1)) print('Negative Eignevalue')
  print(mean(out1[2,]<=0.05))
  file.name<-paste('Exp_norm_4_lowvar_', nl,'.rds',sep='')
  saveRDS(out1, file=file.name)
}

### simulation part 4 Normal short big var ################
for(nl in c(20,50,100,200))
{
  print(paste ('Sample size', nl))
  rep.N<-10^4;set.seed(20250808);
  mu.val.x<-mu.val.y<-c(9,7,4,6); ml<-2*nl;  k<-length(mu.val.x);
  sd.val.x<-sd.val.y<-rep(10,k)
  out1<-simu.norm(rep.N, nl=nl);
  if(any(out1[2,]>1)) print('Negative Eignevalue')
  print(mean(out1[2,]<=0.05))
  file.name<-paste('Exp_norm_4_bigvar_', nl,'.rds',sep='')
  saveRDS(out1, file=file.name)
}


######### post simulation summary ########
## the same part ##
## @knitr read.simu.1
path<-"./Aug_8_Result"
fn.temp<-list.files(path = path);

temp1<-sapply(fn.temp, function(x) strsplit(x,split='.', fixed=TRUE)[[1]][1])
temp2.1<-sapply(temp1[1:10], function(x) strsplit(x,split='_', fixed=TRUE)[[1]][2:5])
temp2.2<-sapply(temp1[11:20], function(x) strsplit(x,split='_', fixed=TRUE)[[1]][2:4])


fn<-paste(path, fn.temp, sep='/');
temp3<-NULL;
for(file in fn)
{
  dat.temp<-readRDS(file)
  temp3<-rbind(temp3, c(mean(dat.temp[2,]<=0.01),mean(dat.temp[2,]<=0.05),mean(dat.temp[2,]<=0.1)))
  print(mean(dat.temp[2,]))
  }
temp.ord<-c(2,4,1,3,5)
temp4.norm<-cbind(t(temp2.1), temp3[1:10,])[c(temp.ord, temp.ord+5),][,-(1:3)]
temp4.exp<-cbind(t(temp2.2), temp3[11:20,])[c(temp.ord, temp.ord+5),][,-(1:2)]
colnames(temp4.norm)<-colnames(temp4.exp)<-c("Sub-Sample Size", '<=0.01','<=0.05', '<=0.1')

saveRDS(temp4.norm, file='norm_simu_same.rds' )
saveRDS(temp4.exp, file='exp_simu_same.rds' )


## the Diff part ##
## @knitr read.simu.2
path<-"./Aug_12_Result"
fn.temp<-list.files(path = path);

temp1<-sapply(fn.temp, function(x) strsplit(x,split='.', fixed=TRUE)[[1]][1])
temp2.1<-sapply(temp1[1:8], function(x) strsplit(x,split='_', fixed=TRUE)[[1]][2:4])
temp2.2<-sapply(temp1[9:20], function(x) strsplit(x,split='_', fixed=TRUE)[[1]][2:4])


fn<-paste(path, fn.temp, sep='/');
temp3<-NULL;
for(file in fn)
{
  dat.temp<-readRDS(file)
  temp3<-rbind(temp3, c(mean(dat.temp[2,]<=0.01),mean(dat.temp[2,]<=0.05),mean(dat.temp[2,]<=0.1)))
}
temp.ord<-c(2,4,1,3)
temp4.exp<-cbind(t(temp2.1), temp3[1:8,])[c(temp.ord, temp.ord+4),][,-(1:2)]
temp4.norm<-cbind(t(temp2.2),temp3[9:20,])[c(temp.ord, temp.ord+4,temp.ord+8),][,-(1:2)]
colnames(temp4.norm)<-colnames(temp4.exp)<-c("Sub-Sample Size", '<=0.01','<=0.05', '<=0.1')

saveRDS(temp4.norm, file='norm_diff.rds' )
saveRDS(temp4.exp, file='exp_diff.rds' )




