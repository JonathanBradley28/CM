#install.packages("devtools") #uncomment this line if you have not installed "devtools"

library(devtools)

#the author token will be removed once the repository  goes public.

install_github("JonathanBradley28/CM",auth_token="7fcfcc652f9d3b1441a8050c13d44f50abb65eaa")

library(CM)

#The help file contains simulation examples.

#To access these examples, run the following commands

help(GibbsBernoulliMLB)

help(GibbsBinomialMLB)

help(GibbsMultinomialMLB)

help(GibbsPoissonMLB)
