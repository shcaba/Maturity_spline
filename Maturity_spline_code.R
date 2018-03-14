library(reshape2)
library(ggplot2)
library(TTR)
library(mgcv)

#Maturity function found in Stock Synthesis
maturity.fxn.SS<-function(lts,slope,Lmat50)
{
  maturity=1/(1+exp(slope*(lts-Lmat50)))
  return(cbind(lts,maturity))
}

#Mat.bins.props
#This function calculates the proportion maturity given
#1) Lengths or ages (data.in)
#2) Number of desired bins (num.bins)
#3) Funtional or biological maturity (fxn_or_bio; FUNCTIONAL=1; BIOLOGICAL=2)
#Returns proportion maturity by length/age bins
Mat.bins.props<-function(data.in,num.bins,fxn_or_bio=1)
{
  name.temp<-names(data.in)[1]
  data.in<-na.omit(data.in)
  names(data.in)[1]<-"METRIC"
  bins<-pretty(data.in$METRIC,num.bins)
  bins_ind<-.bincode(data.in$METRIC,bins)
  data.in$BINS<-bins[bins_ind]
  if(fxn_or_bio==1){Prop_mat_bins<-dcast(data.in,1~METRIC,sum,na.rm=T,value.var = "Functional_maturity")[-1]/table(data.in$METRIC)}
  if(fxn_or_bio==2){Prop_mat_bins<-dcast(data.in,1~METRIC,sum,na.rm=T,value.var = "Biological_maturity")[-1]/table(data.in$METRIC)}
  Data.in.props<-as.data.frame(cbind(as.numeric(colnames(Prop_mat_bins)),t(Prop_mat_bins),table(data.in$METRIC)))
  colnames(Data.in.props)[2:3]<-c("Prop_mat","N")
  print(ggplot(Data.in.props,aes(V1,Prop_mat))+geom_point(aes(size=N))+xlab(paste0(name.temp,"_bins"))) #basic plot to see bin structure#
  colnames(Data.in.props)[1]<-paste0(name.temp,"_bins")
  return(Data.in.props)
}


knot.test.props<-function(knots.in,dat.in)
{
  
  name.in<-names(dat.in)[1]
  dat.in<-na.omit(dat.in)
  names(dat.in)[1]<-"METRIC"
  
  Mat50.out<-CV.crit<-rep(NA,length(knots.in))
  bin.prop.out<-matrix(NA,length(knots.in),length(dat.in$METRIC))
  for(i in 1:length(knots.in))
  {
    spline.out.temp<-smooth.spline(dat.in$METRIC,dat.in$Prop_mat,dat.in$N,all.knots = FALSE,nknots = knots.in[i],cv=TRUE)
    Mat50.out[i]<-uniroot(function(xx) predict(spline.out.temp,xx, type="response")$y - 0.5,range(dat.in$METRIC))$root
    bin.prop.out[i,]<-spline.out.temp$y
    CV.crit[i]<-spline.out.temp$cv.crit
  }
  knots.Mat50.out<-as.data.frame(cbind(knots.in,CV.crit,Mat50.out))
  colnames(knots.Mat50.out)<-c("knots","CrossVal",paste0(name.in,"_50%"))
 #Diagnostic plot
  #Cross validation plot
  print(ggplot(knots.Mat50.out,aes(knots,CrossVal))+geom_point(color="darkgreen",size=4)+labs(x="Knots",y="Cross validation score",title="Cross validation plot"))
  #Length at 50% maturity
  print(ggplot(knots.Mat50.out,aes(knots,get(names(knots.Mat50.out)[3])))+geom_point(color="purple",size=4)+labs(x="Knots",y="Maturity at 50%",title="Maturity at 50% plot")) #L50% plot
  #Proportions at bin plot
  curves.out<-bin.prop.out      #knots.out[,4:ncol(knots.out)]
  rownames(curves.out)<-knots.Mat50.out[,1]
  curves.out.gg<-melt(t(curves.out))
  colnames(curves.out.gg)<-c("bins","knots","props")
  print(ggplot(curves.out.gg,aes(knots,props,color=bins,group=bins))+geom_line()+geom_point()+ 
    ggtitle("Proportions by bins by knots") + xlab("Knots")+ylab("Proportion mature"))
  return(cbind(knots.Mat50.out,bin.prop.out))
}

#Test knots using raw maturity data
knot.test.binary<-function(knots.in,dat.in)
{
  
  name.in<-names(dat.in)[1]
  dat.in<-na.omit(dat.in)
  names(dat.in)[1]<-"METRIC"
  
  Mat50.out<-CV.crit<-rep(NA,length(knots.in))
  for(i in 1:length(knots.in))
  {
    spline.out.temp<-smooth.spline(dat.in$METRIC,dat.in$Functional_maturity,all.knots = FALSE,nknots = knots.in[i],cv=TRUE)
    Mat50.out[i]<-uniroot(function(xx) predict(spline.out.temp,xx, type="response")$y - 0.5,range(dat.in$METRIC))$root
    CV.crit[i]<-spline.out.temp$cv.crit
  }
  knots.Mat50.out<-as.data.frame(cbind(knots.in,CV.crit,Mat50.out))
  colnames(knots.Mat50.out)<-c("knots","CrossVal",paste0(name.in,"_50%"))
  #Cross validation plot
  print(ggplot(knots.Mat50.out,aes(knots,CrossVal))+geom_point(color="blue",size=4)+labs(x="Knots",y="Cross validation score",title="Cross validation plot"))
  #Length at 50% maturity
  print(ggplot(knots.Mat50.out,aes(knots,get(names(knots.Mat50.out)[3])))+geom_point(color="orange",size=4)+labs(x="Knots",y="Maturity at 50%",title="Maturity at 50% plot")) #L50% plot
  return(knots.Mat50.out)
}


Maturity.comp.plots<-function(mat.dat.in,matvals.glm,spline.model)
{
  #mat.nls<-as.data.frame(maturity.fxn.SS(mat.dat.in$Lt_bins,matvals.nls[1],matvals.nls[2]))
  mat.glm<-as.data.frame(maturity.fxn.SS(mat.dat.in$Lt_bins,matvals.glm[1],matvals.glm[2]))
  mat.spline<-as.data.frame(predict(spline.model,mat.dat.in$Lt_bins,type="response"))
  mat.spline$y[mat.spline$y<0]<-0
  mat.spline$y[mat.spline$y>1]<-1
  
  colnames(mat.spline)<-c("lts","maturity")
  
 # mat.nls$model<-"nls"
  mat.glm$model<-"glm"
  mat.spline$model<-"spline"
  
  #mat.nls$N<-mat.glm$N<-mat.spline$N<-YE.mat.props$N
  mat.out<-rbind(mat.glm,mat.spline)
  mat.out$model<-as.factor(mat.out$model)
  mat.gg<-ggplot(mat.out,aes(lts,maturity,color=model))+geom_line(lwd=2)+geom_point(aes(Lt_bins,Prop_mat,size=N),mat.dat.in,color="black")
  print(mat.gg)
  return(mat.spline)
}


################################
################################

######################
### Aurora example ###
######################

#Load data
setwd("D:/JMC/Documents/GitHub/Maturity_spline/")
load("AU_cert.DMP")

#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
Data.in.age<-Data.in[,c(13,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,24,fxn_or_bio=1)
Data.in.bins.age<-Mat.bins.props(Data.in.age,60,fxn_or_bio=1)
#Choose knots for spline; includes cross validation plot as standard output
knots.out<-knot.test.props(knots.in=c(5:24),dat.in=Data.in.bins.lt)
knots.out.fxn<-knot.test.binary(knots.in=c(5:24),dat.in=Data.in.lt)


### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.props,matvals.glm,spline.out.13) ###Adjusting the bins, gives a false results for  N now###






###glm fit###
logistic.mat.fit<-function(Data.in.lt)
fit.mat.glm <- glm (maturity ~ 1 + length, data <-data.frame(length = Data.in.lt$Length, maturity <- Data.in.bins.lt$Functional_maturity),
                    family = binomial(link ="logit"))
summary(fit.mat.glm)
Aglm<- -13.40066
Bglm<- 0.56006
matvals.glm<-c(-fit.mat.glm$coefficients[2], fit.mat.glm$coefficients[1]/-fit.mat.glm$coefficients[2])


#Choose knots
spline.out.13<-smooth.spline(x=Data.in.bins.lt$Length_bins,y=Data.in.bins.lt$Prop_mat,w=Data.in.bins.lt$N,all.knots = FALSE,nknots = 13)
spline.out.9<-smooth.spline(x=Data.in.bins.lt$Length_bins,y=Data.in.bins.lt$Prop_mat,w=Data.in.bins.lt$N,all.knots = FALSE,nknots = 9)
uniroot(function(xx) predict(spline.out.13,xx, type="response")$y - 0.5,range(Data.in.bins.lt$Length_bins))$root ####L50 result###
spline.out.13.01<-smooth.spline(x=Data.in.lt$Length,y=Data.in.lt$Functional_maturity,all.knots = FALSE,nknots = 13)
spline.out.9.01<-smooth.spline(x=Data.in.lt$Length,y=Data.in.lt$Functional_maturity,all.knots = FALSE,nknots = 9)
uniroot(function(xx) predict(spline.out.13.01,xx, type="response")$y - 0.5,range(Data.in.lt$Length))$root ####L50 result###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-as.data.frame(predict(spline.out.13,st.ass.bins,type="response"))
mat.spline.fxn<-as.data.frame(predict(spline.out.13.01,st.ass.bins,type="response"))
plot(mat.spline)
lines(mat.spline.fxn)

mat.spline<-as.data.frame(predict(spline.out.9,st.ass.bins,type="response"))
mat.spline.fxn<-as.data.frame(predict(spline.out.9.01,st.ass.bins,type="response"))
plot(mat.spline)
lines(mat.spline.fxn)



              #####Biological maturity all years -Glm model - calculate 95% CI########
