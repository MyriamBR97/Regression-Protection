# Packages/Libraries installed--------------------------------------------
install.packages("corrplot")
library(corrplot)

install.packages("Rcpp")

devtools::install_github("kassambara/factoextra",force=TRUE)
if(!require(devtools)) install.packages("devtools")

if(!require(ggplot2)) install.packages("ggplot2")
library("ggplot2")

# PREPARATION DE LA BASE DE DONNEES:--------------------------------------

df = read.csv("diametre_asteroid.csv")  # Notre base de donnees
head(df) 

Data<-subset(df,select=-c(full_name,IR,n_obs_used,ad))  #On supprime les colonnes que nous n'allons pas utiliser (full_name,IR,,n_obs_used,ad)
summary(Data)

which (is.na(Data))   #On verifie si il reste des valeurs NULL/NA dans notre base et on supprime les lignes concernees
new_data<-na.omit(Data)
head(new_data)

# REGRESSION MULTIPLE:-----------------------------------------------------

reg.mul<-lm(diameter~.,new_data)  #On applique une regression multiple sur notre base avec diameter la variable a expliquer
summary(reg.mul)
#significative: H, BV puis data_arc

confint(reg.mul)              #intervalles de confiance
reg.mul$residuals
hist(reg.mul$residuals)
shapiro.test(reg.mul$residuals)         #  les variables ne sont pas Gaussiennes car p-value < 5%

reg.mulsansC<-lm(diameter~.-1,new_data)
summary(reg.mulsansC)
reg.mulsansC$residuals
hist(reg.mulsansC$residuals)
shapiro.test(reg.mulsansC$residuals)

reg.mulfin<-lm(diameter~e+i+q+per_y+data_arc+H+albedo+BV+UB+moid-1,new_data)
summary(reg.mulfin)

res.m<-rstudent(reg.mulfin)
plot(res.m,pch=15,cex=.5,ylab="Residus",main="",ylim=c(-3,+3))
abline(h=c(-2,0,2),lty=c(2,1,2))   # On a bien 95% des residus dans la marge
# Le model est suppose correct
mcor<-cor(new_data) #matrice de correlation
mcor

corrplot(mcor,type="upper",order="hclust",tl.col="black",tl.srt=45)  #On remarque une forte correlation entre la variable H et le diametre

xnew<-matrix(c(3,0.14,7.6,160,177,2.4,4.55,46288,8.8,0.094,11,0.72,0.36,1.37),nrow=1)
colnames(xnew)<-c("a","e","i","om","w","q","per_y","data_arc","H","albedo","rot_per","BV","UB","moid") 
xnew<-as.data.frame(xnew)
predict(reg.mulfin,xnew,interval="pred")   # prediction = 140.86 pour intervalle Large [56.04;225.68] car petit echantillon 
#REGRESSION fORWARD:--------------------------------------------
#Procedure manuelle
lm.full<-lm(diameter~.,new_data)
lm.full

lm.null<-lm(diameter~1,new_data)
lm.null

add1(lm.null,scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F")
#la plus petite valeur, donc significative
#c'est la variable a introduire dans le modele.
#ici c'est H
add1(update(lm.null,~.+H),scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F")
#ici c'est BV et albedo
#de faÃ§on totalement aleatoire ,on choisi albedo
add1(update(lm.null,~.+H+albedo),scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F")
#ici c'est data_arc
add1(update(lm.null,~.+H+albedo+data_arc),scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F")
#BV
add1(update(lm.null,~.+H+albedo+data_arc+BV),scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F")
#UB
add1(update(lm.null,~.+H+albedo+data_arc+BV+UB),scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F")
#w
add1(update(lm.null,~.+H+albedo+data_arc+BV+UB+w),scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F")
#a
add1(update(lm.null,~.+H+albedo+data_arc+BV+UB+w+a),scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F")
#per_y
add1(update(lm.null,~.+H+albedo+data_arc+BV+UB+w+a+per_y),scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F")
reg.mul<-lm(diameter~H+albedo+data_arc+BV+UB+w+a+per_y,new_data)
summary(reg.mul)



#forward automatique...juste pour verifier

full<-lm(diameter~.,new_data)
null<-lm(diameter~1,new_data)
forw<-step(null,scope=list(lower=null,upper=full),direction="forward",trace=0)
formula(forw)
summary(forw)

# REGRESSION BACKWARD:------------------------------------------
# Procedure manuelle

lm.full<-lm(diameter~.,new_data) 
drop1(lm.full,test="F")
# on retire du model, les variables les moins significatives une par une
drop1(update(lm.full,~.-moid),test="F")
drop1(update(lm.full,~.-moid-om),test="F") 
drop1(update(lm.full,~.-moid-om-e),test="F")
drop1(update(lm.full,~.-moid-om-e-q),test="F")
drop1(update(lm.full,~.-moid-om-e-q-i),test="F")
drop1(update(lm.full,~.-moid-om-e-q-i-rot_per),test="F")
drop1(update(lm.full,~.-moid-om-e-q-i-rot_per-w),test="F")   # il ne reste plus que des variables significatives
summary(update(lm.full,~.-moid-om-e-q-i-rot_per-w))          # R2 = 0.7334 et p-value < 5%

#Procedure automatique
back<-step(lm.full,direction="backward")
summary(back)                                         # R2 = 0.7351 et p-value < 5%


#REGRESSION SUR ACP:------------------------------------------

cardata <- write.csv(new_data,"acp_data.csv",row.names=FALSE)
cardata <- read.csv("acp_data.csv")
head(cardata)

install.packages("FactoMineR")
library("FactoMineR")
library("factoextra")

res.pca=PCA(cardata[,1:14],scale.unit=TRUE,ncp=14,graph=T)
#Les variables positivement corrélées sont regroupées.
#Les variables négativement corrélées sont positionnées sur les côtés opposés de l'origine du graphique.
#La distance entre les variables et l'origine mesure la qualité de représentation des variables. 
#Les variables qui sont loin de l'origine sont bien représentées par l'ACP.

var <- get_pca_var(res.pca)
head(var$cos2)      # Resultats de l'ACP

fviz_pca_var(res.pca, col.var = "cos2",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE) # Coloration et pas de chevauchement de texte


res.pca$eig   #valeur propre de la matrice de correlation

fviz_eig(res.pca, geom="bar", width=0.8, addlabels=T)

res.pca$ind$coord
novariables<-res.pca$ind$coord    #matrice des scores
novariables                       #vecteurs othogonaux


#REGRESSION FORWARD POUR ACP:------------------------------------------
#Procedure manuelle

full<-lm(diameter~.,data=cardata)
full   # beta0= 1.297e+03, beta1= -4.581e+02,.... 

null<-lm(diameter~1,data=cardata)

add1(null,scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F") #On cherche la variable la + correlee au diametre donc la F-Value la plus grande = H
add1(update(null,~.+H),scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F") # donc la F-Value la plus grande = BV
reg.mul_ACP<-lm(diameter~H,data=cardata)
summary(reg.mul_ACP)

add1(update(null,~.+H+BV),scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F") #donc la F-Value la plus grande = data_arc
add1(update(null,~.+H+BV+data_arc),scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F") #donc la F-Value la plus grande = UB
add1(update(null,~.+H+BV+data_arc+UB),scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F") #donc la F-Value la plus grande = albedo
add1(update(null,~.+H+BV+data_arc+UB+albedo),scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F") #donc la F-Value la plus grande = a
add1(update(null,~.+H+BV+data_arc+UB+albedo+a),scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F") #donc la F-Value la plus grande = per_y
add1(update(null,~.+H+BV+data_arc+UB+albedo+a+per_y),scope=~a+e+i+om+w+q+per_y+data_arc+H+albedo+rot_per+BV+UB+moid,test="F") #il ne reste plus de variables significatives

reg.mul_ACP<-lm(diameter~H+BV+data_arc+UB+albedo+a+per_y,data=cardata)
summary(reg.mul_ACP)

#Procedure automatique
forw<-step(null,scope=list(lower=null,upper=full),direction="forward",trace=0)
formula(forw)
summary(lm(forw,data=cardata))

#REGRESSION RIDGE:------------------------------------------
library(MASS)
mod.ridge=lm.ridge(diameter~.,new_data,lambda=seq(0,4,0.01))
plot(mod.ridge)
legend("right",legend=rownames(mod.ridge$coef),cex=0.5,col=1:14,lty=1:3)

select(mod.ridge)
mod.ridgeopt=lm.ridge(diameter~.,new_data,lambda=0.12)
coeff=coef(mod.ridgeopt)
coeff

#_____________FIN_______________
