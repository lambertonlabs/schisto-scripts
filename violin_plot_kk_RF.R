dat_BW_KK = read.csv("BW_KK_summary.csv", header= T, stringsAsFactors = FALSE)



ggplot(dat_BW_KK,aes(x=as.factor(timept),y=mean)) +geom_violin() +scale_y_log10()

ggplot(dat_BW_KK,aes(x=as.factor(timept),y=mean)) +geom_boxplot() +scale_y_log10() 
  + scale_fill_discrete(name="TimePoint", labels=c("baseline", "5months", "3wks postTx"))
