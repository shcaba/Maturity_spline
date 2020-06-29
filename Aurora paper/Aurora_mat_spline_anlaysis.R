######################
### Aurora example ###
######################

#Load data
setwd("//nwcfile/FRAM/Users/Melissa.Head/Aurora_Mat/Aurora manuscript")
AU.mat<-read.csv("2012_2014_WCGBT_AuroramaturityEdited.csv")
AU.cert<-subset(AU.mat,Certainty==1)
save(AU.mat,file="AU_mat.DMP")
AU.age <-read.csv("Aurora2012_propage.csv")
AU.age.cert<-subset(AU.age,Certainty==1)
load("AU_mat.DMP")


###########Biological length at maturity - all years#######################
Data.in<-AU.cert


#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=2)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:24),dat.in=Data.in.lt) #15+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:24),dat.in=Data.in.bins.lt) #13+ looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Biological_maturity,all.knots = FALSE,nknots = 15)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,2,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")

#################Biological maturity - Pass 1#########
Data.in<-AU_cert[AU_cert$Pass==1,]


#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=2)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:24),dat.in=Data.in.lt) #15+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:24),dat.in=Data.in.bins.lt) #13+ looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Biological_maturity,all.knots = FALSE,nknots = 15)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,2,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")

#################Functional maturity - Pass 1#########
Data.in<-AU.cert[AU.cert$Pass==1,]


#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=1)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:24),dat.in=Data.in.lt) #16+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:24),dat.in=Data.in.bins.lt) #13+ looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model

spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Functional_maturity,all.knots = FALSE,nknots = 13)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,3,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)



#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 13)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")



#################Biological maturity - Pass 2#########
Data.in<-AU_cert[AU_cert$Pass==2,]


#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=2)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:24),dat.in=Data.in.lt) #15+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:24),dat.in=Data.in.bins.lt) #13+ looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Biological_maturity,all.knots = FALSE,nknots = 15)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,2,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")

#################Functional maturity - Pass 2#########
Data.in<-AU.cert[AU.cert$Pass==2,]


#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=1)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:24),dat.in=Data.in.lt) #16+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:24),dat.in=Data.in.bins.lt) #13+ looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model

spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Functional_maturity,all.knots = FALSE,nknots = 13)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,3,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)



#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 13)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")




###########Functional length at maturity - all years#######################
Data.in<-AU.cert


#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=1)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:24),dat.in=Data.in.lt) #15+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:24),dat.in=Data.in.bins.lt) #16 looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Functional_maturity,all.knots = FALSE,nknots = 16)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,3,1000,11,Lmat_50)
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")


###########Biological length at maturity - 2012#######################

Data.in<-subset(AU.cert,Year==2012)


#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=2)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:20),dat.in=Data.in.lt) #9+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:20),dat.in=Data.in.bins.lt) #11+ looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Biological_maturity,all.knots = FALSE,nknots = 11)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,2,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")


###########Functional length at maturity - 2012#######################

Data.in<-subset(AU.cert,Year==2012)

#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=1)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:20),dat.in=Data.in.lt) #10+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:20),dat.in=Data.in.bins.lt) #8+ looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Functional_maturity,all.knots = FALSE,nknots = 10)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,3,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")


###########Biological length at maturity - 2013#######################

Data.in<-subset(AU.cert,Year==2013)


#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=2)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:20),dat.in=Data.in.lt) #7+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:20),dat.in=Data.in.bins.lt) #8+ looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Biological_maturity,all.knots = FALSE,nknots = 8)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,2,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")


###########Functional length at maturity - 2013#######################

Data.in<-subset(AU.cert,Year==2013)

#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=1)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:20),dat.in=Data.in.lt) #7+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:20),dat.in=Data.in.bins.lt) #8+ looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Functional_maturity,all.knots = FALSE,nknots = 8)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,3,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")

###########Biological length at maturity - 2014#######################

Data.in<-subset(AU.cert,Year==2014)


#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=2)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:20),dat.in=Data.in.lt) #11+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:20),dat.in=Data.in.bins.lt) #8+ looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Biological_maturity,all.knots = FALSE,nknots = 11)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,2,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")


###########Functional length at maturity - 2014#######################

Data.in<-subset(AU.cert,Year==2014)

#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=1)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:20),dat.in=Data.in.lt) #11+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:20),dat.in=Data.in.bins.lt) #10+ looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Functional_maturity,all.knots = FALSE,nknots = 11)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,3,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,3,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")

###########Biological length at maturity - 2015#######################

Data.in<-subset(AU.cert,Year==2015)


#Prep data: age or length
Data.in.lt<-Data.in[,c(9,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=2)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:15),dat.in=Data.in.lt) #8-10 looks good
knots.out.prop<-knot.test.props(knots.in=c(5:20),dat.in=Data.in.bins.lt) #8+ looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Biological_maturity,all.knots = FALSE,nknots = 9)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,2,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")


###########Functional length at maturity - 2015#######################

Data.in<-subset(AU.cert,Year==2015)

#Prep data: age or length
Data.in.lt<-Data.in[,c(10,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=1)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:20),dat.in=Data.in.lt) #5, 10-14 looks good
knots.out.prop<-knot.test.props(knots.in=c(5:20),dat.in=Data.in.bins.lt) #5, 9-14 looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Functional_maturity,all.knots = FALSE,nknots = 16)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,3,10000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,3,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")


###########Biological length at maturity - 2016#######################

Data.in<-subset(AU.cert,Year==2016)


#Prep data: age or length
Data.in.lt<-Data.in[,c(10,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=2)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:20),dat.in=Data.in.lt) #11+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:20),dat.in=Data.in.bins.lt) #10+ looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Biological_maturity,all.knots = FALSE,nknots = 11)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,2,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")


###########Functional length at maturity - 2016#######################

Data.in<-subset(AU.cert,Year==2016)

#Prep data: age or length
Data.in.lt<-Data.in[,c(10,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=1)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:20),dat.in=Data.in.lt) #5, 11+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:20),dat.in=Data.in.bins.lt) #9,15,  looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Functional_maturity,all.knots = FALSE,nknots = 15)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,3,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")

####################North of Mendocino - binned above 37 cm, Biological maturity#####



Data.in <- AU.cert[AU.cert$Latitude <= 48.43 & AU.cert$Latitude > 40.433,]


#Prep data: age or length
Data.in.lt<-Data.in[,c(10,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=2)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:20),dat.in=Data.in.lt) #10-15 looks good
knots.out.prop<-knot.test.props(knots.in=c(5:20),dat.in=Data.in.bins.lt) #10-15 looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Biological_maturity,all.knots = FALSE,nknots = 11)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,2,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")


###########North of Mendocino Functional length at maturity, binned above 37cm, all years#######################

Data.in <- AU.cert[AU.cert$Latitude <= 48.43 & AU.cert$Latitude > 40.433,]

#Prep data: age or length
Data.in.lt<-Data.in[,c(10,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=1)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:20),dat.in=Data.in.lt) # 10-15 looks good
knots.out.prop<-knot.test.props(knots.in=c(5:20),dat.in=Data.in.bins.lt) #10-15  looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Functional_maturity,all.knots = FALSE,nknots = 10)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,3,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,3,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")


####################South of Mendocino - binned above 37 cm, Biological maturity#####


Data.in <- AU.cert[AU.cert$Latitude <=40.433 & AU.cert$Latitude > 32.0,]

#Prep data: age or length
Data.in.lt<-Data.in[,c(10,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=2)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:20),dat.in=Data.in.lt) #10+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:20),dat.in=Data.in.bins.lt) #10+

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Biological_maturity,all.knots = FALSE,nknots = 11)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,2,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")


###########South of Mendocino Functional length at maturity, binned above 37cm, all years#######################

Data.in <- AU.cert[AU.cert$Latitude <=40.433 & AU.cert$Latitude > 32.0,]

#Prep data: age or length
Data.in.lt<-Data.in[,c(10,15,16)]
#Define bin structure
Data.in.bins.lt<-Mat.bins.props(Data.in.lt,seq(0,40,1),fxn_or_bio=1)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi<-knot.test.binary(knots.in=c(5:20),dat.in=Data.in.lt) # 11+ looks good
knots.out.prop<-knot.test.props(knots.in=c(5:20),dat.in=Data.in.bins.lt) #13+  looks good

#Fit binary data models
Data.in.model<-Data.in.lt
###glm fit###
glm.model<-logistic.mat.fit(Data.in.model)
#Spline model
spline.out<-smooth.spline(x=Data.in.model$Length,y=Data.in.model$Functional_maturity,all.knots = FALSE,nknots = 10)
Lmat_50<-uniroot(function(xx) predict(spline.out,xx, type="response")$y - 0.5,range(Data.in.model$Length))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_boot<-Boot.Lmat50(Data.in.model,1,3,1000,11,Lmat_50)####change between functional and biological data input###
quantile(Lmat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)


#Fit spline model using proportions
Data.in.props.model<-Data.in.bins.lt
spline.out.bins<-smooth.spline(x=Data.in.props.model$Length_bins,y=Data.in.props.model$Prop_mat,w=Data.in.props.model$N,all.knots = FALSE,nknots = 16)
Lmat_50_bins<-uniroot(function(xx) predict(spline.out.bins,xx, type="response")$y - 0.5,range(Data.in.props.model$Length_bins))$root ####L50 result###
#Calculate uncertainty in Lmat50%
Lmat_50_bins_boot<-Boot.Lmat50(Data.in.props.model,1,2,1000,11,Lmat_50_bins)
quantile(Lmat_50_bins_boot,probs=c(0.025,0.5,0.975),na.rm=TRUE)

### Plot comparisons ###
matplot<-Maturity.comp.plots(Data.in.model,Data.in.props.model,seq(0,38,1),glm.model$parameters,spline.out,spline.out.bins) ###Adjusting the bins, gives a false results for  N now###

#Predict spline for assessment
st.ass.bins<-seq(8,40,2)
mat.spline<-Spline.fit(spline.out,st.ass.bins,data.type="Lengths")
plot(mat.spline,type="b")



###################
### Age example ###
###################

###############Age at maturity 2012 only -  biological maturity##############
 
Data.in <-AU.age.cert
Data.in.age<-Data.in[,c(11,13,14)]
Data.in.age<-na.omit(Data.in.age)
Data.in.bins.age<-Mat.bins.props(Data.in.age,c(1:80),fxn_or_bio=2)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi.age<-knot.test.binary(knots.in=c(5:20),dat.in=Data.in.age) #10+ looks good
knots.out.prop.age<-knot.test.props(knots.in=c(5:24),dat.in=Data.in.bins.age) #13+ looks good

#Fit binary data models
Data.in.model.ages<-Data.in.age
###glm fit###
glm.model.ages<-logistic.mat.fit(Data.in.model.ages)
#Spline model
spline.out.ages<-smooth.spline(x=Data.in.model.ages$Age,y=Data.in.model.ages$Biological_maturity,all.knots = FALSE,nknots = 14)
Amat_50<-uniroot(function(xx) predict(spline.out.ages,xx, type="response")$y - 0.5,range(Data.in.model.ages$Age))$root ####L50 result###

Amat_50_boot<-Boot.Lmat50(Data.in.age,1,2,10000,11,Amat_50) ###make
quantile(Amat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)
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

###############Age at maturity 2012 -  functional maturity####
 
Data.in <-AU.age.cert
Data.in.age<-Data.in[,c(11,13,14)]
Data.in.age<-na.omit(Data.in.age)
Data.in.bins.age<-Mat.bins.props(Data.in.age,c(1:80),fxn_or_bio=1)
#Choose knots for spline; includes cross validation plot as standard output
knots.out.bi.age<-knot.test.binary(knots.in=c(5:20),dat.in=Data.in.age) #13+ looks good
knots.out.prop.age<-knot.test.props(knots.in=c(5:24),dat.in=Data.in.bins.age) #13+ looks good

#Fit binary data models
Data.in.model.ages<-Data.in.age
###glm fit###
glm.model.ages<-logistic.mat.fit(Data.in.model.ages)
#Spline model
spline.out.ages<-smooth.spline(x=Data.in.model.ages$Age,y=Data.in.model.ages$Functional_maturity,all.knots = FALSE,nknots = 14)
Amat_50<-uniroot(function(xx) predict(spline.out.ages,xx, type="response")$y - 0.5,range(Data.in.model.ages$Age))$root ####L50 result###

Amat_50_boot<-Boot.Lmat50(Data.in.age,1,3,10000,11,Amat_50) ###make
quantile(Amat_50_boot,probs=c(0.05,0.5,0.95),na.rm=TRUE)
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





