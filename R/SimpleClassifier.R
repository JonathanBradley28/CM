SimpleClassifier<-function(data,phat,holdout=NULL,phathold=NULL){

  cutoff<-function(lambda){
    estpos=phat>lambda
    estneg = phat<=lambda
    truepos=data==1
    trueneg=data==0

    falsepositive= mean(estpos[trueneg])
    falsenegative= mean(estneg[truepos])

    helliginer= sqrt(falsepositive)+sqrt(falsenegative)
        return(helliginer)
  }

  initB=0.5
  estpar<-optim(initB,cutoff)
  estpos.final=phat>estpar$par[1]

  helliginer=NULL
  estpos.hold=NULL
if(length(holdout)>0){
  truepos=holdout==1
  trueneg=holdout==0

  estpos.hold=phathold>estpar$par[1]
  estneg.hold=phathold<=estpar$par[1]

  falsepositive= mean(estpos.hold[trueneg])
  falsenegative= mean(estneg.hold[truepos])

  helliginer= sqrt(falsepositive)+sqrt(falsenegative)
}

  output<-list(helliginer,estpar$par[1],estpos.final,estpos.hold)

  return(output)
}
