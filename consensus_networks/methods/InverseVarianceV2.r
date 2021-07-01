 ###################################################################################################
    ################################################################################################### 
    #                                                                                                 #
    #               Consensus network: Inverse variance weighted method                        #                 
    #                                                                                                 #
    ################################################################################################### 


### Meta-analysis
### Iterate through each edge and combine statistical confidences using
### meta-analysis technique inverse variance-weighted method

MetaAnalysis = function(x) 
{
 pooled_est <- 0;
 weight_sum<-0;
 
 
 ### Number of networks learnt that an edge exist
 m<-length(x)
 
  for (k in 1:length(x))
        {
            d<-x[k]*m;
            if (d == 0){
                dataset_est<-0;
                var<-0;
                weight<-0;
                }
            else {
                dataset_est = log(x[k])
                var = 1/d;
                weight=1/var;
            pooled_est <- pooled_est + weight*dataset_est;
            weight_sum <- weight_sum + weight;
          
            
        }
    }
 
         if(weight_sum == 0){
            MetaAnalysis<-0
            } else         
           MetaAnalysis <- exp(pooled_est/weight_sum)

  return(MetaAnalysis)

 }
 
  ################################################################################################### 
 
red_wg_p<-DF[,2:(dim(DF)[2])]
 

Meta.score=matrix()

 for (i in 1:dim(red_wg_p)[1])
  {
   Meta.score[i]<-MetaAnalysis(red_wg_p[i,])
  }
   
  Meta.score <-unlist(Meta.score)
  
