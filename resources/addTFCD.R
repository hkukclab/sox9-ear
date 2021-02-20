load("../resources/TF.CD.receptors.human.mouse.RData")

#####################################################
Add_TFCD<-function(gin,species="mouse"){
	gin<-as.character(gin)
	if(species=="mouse"){
		ind1<-which(gin%in%TF.mouse)
		ind2<-which(gin%in%surface.mouse)
	}else{
		ind1<-which(gin%in%TF.human)
		ind2<-which(gin%in%surface.human)
	}
	gin[ind1]<-paste0(gin[ind1],"*")
	gin[ind2]<-paste0(gin[ind2],"#")
	gin
}
