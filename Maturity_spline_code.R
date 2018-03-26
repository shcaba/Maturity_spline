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
#2) Number of desired bins OR pre-specified bin structure (num.bins)
#3) Funtional or biological maturity (fxn_or_bio; FUNCTIONAL=1; BIOLOGICAL=2)
#Returns proportion maturity by length/age bins
Mat.bins.props<-function(data.in,num.bins,fxn_or_bio=1)
{
  name.temp<-names(data.in)[1]
  data.in<-na.omit(data.in)
  names(data.in)[1]<-"METRIC"
  if(length(num.bins)==1){bins<-pretty(data.in$METRIC,num.bins)}
  if(length(num.bins)>1){bins<-num.bins}
  bins_ind<-.bincode(data.in$METRIC,bins)
  data.in$BINS<-bins[bins_ind]
  data.in<-na.omit(data.in)
  if(fxn_or_bio==1){Prop_mat_bins<-dcast(data.in,1~BINS,sum,na.rm=T,value.var = "Functional_maturity")[-1]/table(data.in$BINS)}
  if(fxn_or_bio==2){Prop_mat_bins<-dcast(data.in,1~BINS,sum,na.rm=T,value.var = "Biological_maturity")[-1]/table(data.in$BINS)}
  Data.in.props<-as.data.frame(cbind(as.numeric(colnames(Prop_mat_bins)),t(Prop_mat_bins),table(data.in$BINS)))
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

logistic.mat.fit<-function(Data.in)
{
  fit.mat.glm <- glm (maturity ~ 1 + length, data <-data.frame(length = Data.in$Length, maturity <- Data.in$Functional_maturity),
                      family = binomial(link ="logit"))
  
  matvals.glm<-c(-fit.mat.glm$coefficients[2], fit.mat.glm$coefficients[1]/-fit.mat.glm$coefficients[2])
  return(list(logistic.model=fit.mat.glm,parameters=matvals.glm))
}

Maturity.comp.plots<-function(mat.dat.in,mat.props.in,bins.in,logistic.parms,spline.model,spline.model.bins)
{
  #Logistic predictions
  mat.glm<-as.data.frame(maturity.fxn.SS(bins.in,logistic.parms[1],logistic.parms[2]))
  mat.glm$model<-"glm"
  #Spline binary data predictions
  mat.spline<-as.data.frame(predict(spline.model,bins.in,type="response"))
  mat.spline$y[mat.spline$y<0]<-0
  mat.spline$y[mat.spline$y>1]<-1
  colnames(mat.spline)<-c("lts","maturity")
  mat.spline$model<-"spline"
  mat.out<-rbind(mat.glm,mat.spline)
  #Spline binned predictions
  if(exists("spline.model.bins")==T)
  {
    mat.spline.bins<-as.data.frame(predict(spline.model.bins,bins.in,type="response"))
    mat.spline.bins$y[mat.spline.bins$y<0]<-0
    mat.spline.bins$y[mat.spline.bins$y>1]<-1
    colnames(mat.spline.bins)<-c("lts","maturity")
    mat.spline.bins$model<-"spline_bins"
    mat.out<-rbind(mat.glm,mat.spline,mat.spline.bins)
  }
  #Plot models
  mat.out$model<-as.factor(mat.out$model)
  mat.gg<-ggplot(mat.out,aes(lts,maturity,color=model))+geom_line(lwd=2)+geom_point(aes(Length_bins,Prop_mat,size=N),mat.props.in,color="black")
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
load("ARRA_dat.DMP")

Data.in<-AU.cert
#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
Data.in.age<-Data.in[,c(13,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,24,fxn_or_bio=1)
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,2),fxn_or_bio=1)
Data.in.bins.age<-Mat.bins.props(Data.in.age,60,fxn_or_bio=1)
#Choose knots for spline; includes cross validation plot as standard output
knots.out<-knot.test.props(knots.in=c(5:24),dat.in=Data.in.bins.lt)
knots.out.fxn<-knot.test.binary(knots.in=c(5:24),dat.in=Data.in.lt)

Data.in<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in)
#Spline model
spline.out<-smooth.spline(x=Data.in$Length,y=Data.in$Functional_maturity,all.knots = FALSE,nknots = 10)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in$Length))$root ####L50 result###
Data.in.bins<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.bins$Length_bins,y=Data.in.bins$Prop_mat,w=Data.in.bins$N,all.knots = FALSE,nknots = 10)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.bins$Length_bins))$root ####L50 result###

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in,Data.in.bins.lt,seq(0,38,2),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###




#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-as.data.frame(predict(spline.out,st.ass.bins,type="response"))
plot(mat.spline)
lines(mat.spline.fxn)

mat.spline<-as.data.frame(predict(spline.out.9,st.ass.bins,type="response"))
mat.spline.fxn<-as.data.frame(predict(spline.out.9.01,st.ass.bins,type="response"))
plot(mat.spline)
lines(mat.spline.fxn)
