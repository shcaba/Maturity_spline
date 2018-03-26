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
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in$Functional_maturity,all.knots = FALSE,nknots = 11)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###

#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 13)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,36,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="lengths")
plot(mat.spline,type="b")


#Age example
Data.in.age<-Data.in[,c(13,15,16)]
Data.in.bins.age<-Mat.bins.props(Data.in.age,60,fxn_or_bio=1)