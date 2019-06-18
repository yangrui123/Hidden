# Rscript
library(ggplot2)

args = comandArgs(T)

filename = args[1]
dt = read.csv(filename, header = TRUE, sep='\t')
rbps = dt$rbp
spearmans = dt$scc_lncRNA_seq

name = c()
y_pred = c()
y_true = c()


for(i in c(69:79)) {
  filename = paste('lncRNA/scatter_plots/', rbps[i], '.lncRNA.scatter.plot', sep = '')
  dt = read.csv(filename, header = TRUE, sep='\t')
  rbp = paste(rbps[i], '\n', sep = '')
  title = paste(rbp, 'r=', spearmans[i], sep = '')
  subs = rep(title, length(dt$y_pred))
  name = c(name, subs)
  y_pred = c(y_pred, dt$y_pred)
  y_true = c(y_true, dt$y_true)
}

df = data.frame(y_true=y_true, y_pred=y_pred, name=name)
p = ggplot(data=df, aes(x=y_true, y=y_pred))+geom_point(size=1)+
  labs(x='CLIP-seq signal', y='Predicted signal')+
  theme(axis.text = element_text(size=12, face='bold'), text = element_text(size=12, face='bold'),
        strip.text.x = element_text(size = 12, face='bold'))+
  facet_wrap(.~name, scales = "free", nrow=2)

ggsave(p, filename = paste(path, '/figure5a_results/', rbp, '.png', sep=''))
