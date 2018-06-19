######################
### Aurora example ###
######################

#Load data
setwd("D:/JMC/Documents/GitHub/Maturity_spline/")
load("ARRA_dat.DMP")
Data.in<-AU.cert

#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,37,1),fxn_or_bio=1)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:24),dat.in=Data.in.lt) #11 looks good
knots.out.prop<-knot.test.props(knots.in=c(5:24),dat.in=Data.in.bins.lt) #13 looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Functional_maturity,all.knots = FALSE,nknots = 11)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,3,1000,11,Lmat_50)
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 13)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,36,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")

###################
### Age example ###
###################
Data.in.age<-Data.in[,c(13,15,16)]
Data.in.age<-na.omit(Data.in.age)
Data.in.bins.age<-Mat.bins.props(Data.in.age,c(1:80),fxn_or_bio=1)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi.age<-knot.test.binary(knots.in=c(5:24),dat.in=Data.in.age) #11 looks good
knots.out.prop.age<-knot.test.props(knots.in=c(5:24),dat.in=Data.in.bins.age) #13 looks good

#Fit binary data models
Data.in.model.ages<-Data.in.age
###glm fit###
glm.model.ages<-logistic.mat.fit(Data.in.model.ages)
#Spline model
spline.out.ages<-smooth.spline(x=Data.in.model.ages$Age,y=Data.in.model.ages$Functional_maturity,all.knots = FALSE,nknots = 7)
Amat_50<-uniroot(function(xx) predict(spline.out.ages,xx, type="response")$y - 0.5,range(Data.in.model.ages$Age))$root ####L50 result###

#Fit spline model using proportions
Data.in.props.model.ages<-Data.in.bins.age
spline.out.bins.ages<-smooth.spline(x=Data.in.props.model.ages$Age_bins,y=Data.in.props.model.ages$Prop_mat,w=Data.in.props.model.ages$N,all.knots = FALSE,nknots = 13)
Amat_50_bins<-uniroot(function(xx) predict(spline.out.bins.ages,xx, type="response")$y - 0.5,range(Data.in.props.model.ages$Age_bins))$root ####L50 result###

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model.ages,Data.in.props.model.ages,seq(1,80,1),glm.model.ages$parameters,spline.out.ages,spline.out.bins.ages) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(0,80,1)
mat.spline.ages<-Spline.fit(spline.out.ages,st.ass.bins,data.type="Age")
plot(mat.spline.ages,type="b")
