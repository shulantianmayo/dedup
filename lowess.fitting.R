cat("-- reading arguments from ……\n", sep = "");
cmd_args = commandArgs(trailingOnly=TRUE);

for (arg in cmd_args) {
cat("  ", arg, "\n", sep="")
}



fn1=unlist(strsplit(cmd_args[1], "="))[2]
fn2=unlist(strsplit(cmd_args[2], "="))[2]
fn3=unlist(strsplit(cmd_args[3], "="))[2]


data_df=read.delim(fn1, sep="\t", header=F)
colnames(data_df)=c("PEAK_ID", "IP_DUP_perKb", "IP_noDUP_perKb")


my_fit=lowess(as.numeric(data_df$IP_noDUP_perKb), as.numeric(data_df$IP_DUP_perKb), f=as.numeric(fn3))
my_fit_df=data.frame(cbind(my_fit$x, my_fit$y))


data_df_IP=data_df[,c(1,3,2)]
data_df_IP_order=data_df_IP[order(as.numeric(data_df_IP[,2]), decreasing=F),]


my_fit_df_final=cbind(data_df_IP_order, my_fit_df)
colnames(my_fit_df_final)[4:5]=c("IP_non_DUP_perKb", "IP_DUP_perKb_fitted")

write.table(my_fit_df_final, file=fn2, sep="\t", row.names=F, quote=F)

