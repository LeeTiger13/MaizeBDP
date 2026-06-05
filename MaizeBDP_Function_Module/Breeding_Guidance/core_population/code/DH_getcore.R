library(corehunter)
# 导入数据，为遗传距离4DH冗余度检测.py生成的csv格式的文件
data <- distances(file = 'output.DH.dist.txt')
# 提取核心集，需要把size=这个参数设置成用户可以自己输入
core <- sampleCore(a, size = 0.5)
# 输出核心集的样本ID
id <- as.vector(core[[1]])
id_df <- as.data.frame(id)
write.table(id_df[1], file = paste0('selected_id.txt'), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
