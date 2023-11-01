# Statement ----
# 本代码用来分析京都市生物多样性和社会经济因素之间的关系

# Package ---
library(openxlsx)
library(dplyr)
library(tidyr)
library(psych)
library(ggplot2)
library(vegan)
library(geosphere)
library(leaps)
library(sf)
library(terra)
library(tmap)
library(showtext)
showtext_auto()
# if dir "ProcData" does not exist, create one
if(!file.exists("ProcData")) dir.create("ProcData")

# Function ----
# function: get community data based on individual data
# note: output is wide data of community - rows as species and columns as number of trees
# variable: 
# x: investigation data of each tree
# col.group: column name group based on, without quotation marks
GetComm <- function(x, col.group) {
  x %>% 
    mutate(stem = 1) %>%  # 每棵树的丰度即为1
    select({{col.group}}, stem, species) %>%
    pivot_wider(names_from = species, values_from = stem,
                values_fn = sum, values_fill = 0) %>% 
    return()
}

# function: get biodiversity indexes based on individual data 
# variable: 
# x: investigation data of each tree
# x_comm: community data, you can get it with GetComm()
# col.group: column name group based on, without quotation marks
GetDiv <- function(x, x_comm, col.group) {
  # inside function: summary each attribute 
  funin_attrcalc <- function(coltar, tarvalue) {
    x_sub <- x
    x_sub["tarornot"] <- x_sub[coltar] == tarvalue
    x_sub <- x_sub %>% group_by({{col.group}}) %>% 
      summarise(
        perc = sum(1 * tarornot) / sum(1)) %>%
      ungroup() %>% 
      select({{col.group}}, perc)
    names(x_sub)[2] <- paste0("perc_", tarvalue)
    return(x_sub)
  }
  
  output <- x_comm %>%
    mutate(abundance = rowSums(.[3:ncol(.)]),
           richness = apply(.[2:ncol(.)]>0, 1, sum),
           shannon = diversity(.[2:ncol(.)], index = "shannon"),
           simpson = diversity(.[2:ncol(.)], index = "simpson"),
           evenness = shannon / log(richness)) %>%
    select({{col.group}}, 
           abundance, richness, shannon, simpson, evenness)
  
  return(output)
}

# 函数：可视化数据框各列两两之间的相关关系
# 漏洞：暂时以总richness作为指标
# variable: 
# x: independent variable 
# y: dependent variable
GetCorrplot <- function(x, y) {
  # Change colnames. 
  colnames(x) <- x.name
  colnames(y) <- y.name
  # 计算各列两两之间的相关性
  cor.res <- corr.test(x, y)
  # 作图表示相关性大小和是否显著，如果不显著的话，会以打叉表示
  corrplot::corrplot(corr = cor.res$r, method = "number", p.mat = cor.res$p)
}

# 备选变量：仅考虑上一步皮尔森检测中数据量大于10个的变量
GetRegSubset <- function(x, var_response) {
  leaps <- regsubsets(
    richness ~ pop_prop_65_up + pop_prop_75_up + popf_prop_75_up + 
      residential + multi_family_residential + commercial_industrial, 
    data = x)
  plot(leaps, scale = "adjr2", main = var_response)
}

# function: trun lm() result into data.frame
# parameters: 
# x: result of lm()
LmRes2Df <- function(x) {
  # get coefficients 
  ressum <- summary(x)$coefficients
  # get var name, estimate, and p-value 
  output <- data.frame(
    variable = names(ressum[-1, 1]), 
    estimate = ressum[-1, 1], 
    p = ressum[-1, 4]
  ) %>% 
    mutate(p_mark = case_when(
      p < 0.001 ~ "***", 
      p < 0.01 ~ "**", 
      p < 0.05 ~ "*", 
      TRUE ~ ""
    ))
  
  # eliminate rownames
  rownames(output) <- NULL
  
  output$result <- round(output$estimate, digits = 2)
  output$p <- output$p_mark
  output <- output[c("variable", "result", "p_mark")]
  return(output)
}

# Read data ----
## Constant ----
kIndex <- c("abundance", "richness", "shannon", "simpson", "evenness")

# population and structure factors
kPopFac <- c("pop_tot", "pop_prop_0_14", "pop_prop_65_up", "pop_prop_75_up",
            "popf_tot", "popf_prop_0_14", "popf_prop_65_up", "popf_prop_75_up")

kWard <- c("Ukyo-ku", "Kita-ku", "Sakyo-ku",
           "Kamigyo-ku", "Nakagyo-ku", "Shimogyo-ku", "Higashiyama-ku",
           "Minami-ku")

# bug: 是否要仅保留在各个样方中均比较普遍的土地覆盖类型呢？但是如果那样做了，就只剩下residential和transportation两种了。在后期可以考虑在样方周围一定缓冲区内，计算土地利用或覆盖（但是数据不如现场调查那么精确）比例或景观指数和样方内的多样性指数之间的关系？
kLandCover <- c("residential", "transportation", "temple_shrine",
                "multi_family_residential", "agriculture", 
                "commercial_neighbor", "water_wetland", "park", 
                "cemetery", "vacant", "commercial_industrial", 
                "institutional", "golf_course")

# 作图范围
kBBox <- c(135.62, 34.90, 135.839, 35.079)

## Table data ----
# 读取调查数据
indv.tree <- read.csv("RawData/KUP Plant_data.csv") %>% 
  rename_with(tolower) %>% 
  rename(species = species_lt) %>% 
  tibble()
# 聚合成样方级别的数据
qua.bd <- GetComm(indv.tree, plot_id) %>% 
  GetDiv(x = indv.tree, x_comm = ., col.group = plot_id)

# land cover proportion data
land.cover <- read.xlsx("RawData/GIS Quadrat_land_cover.xlsx") %>% 
  as_tibble() %>% 
  rename_with(tolower) %>% 
  rename(qua_id = quadrat_id, land_cover = detailed_land_cover) %>% 
  select(qua_id, land_cover, shape_area) %>% 
  # replace the name with special marks 
  mutate(land_cover = gsub("/", "_", .$land_cover)) %>% 
  mutate(land_cover = gsub(" ", "_", .$land_cover)) %>% 
  mutate(land_cover = gsub("-", "_", .$land_cover)) %>% 
  group_by(qua_id, land_cover) %>% 
  summarise(area = sum(shape_area)/400) %>% 
  pivot_wider(
    id_cols = qua_id, names_from = land_cover, values_from = area) %>% 
  # replace all NA with 0
  mutate(across(everything(), ~replace_na(.x, 0)))

# 社会经济因子
# bug: 有些样地的社会经济因子为NA？
socio.economic <- read.xlsx("RawData/GIS Socio_economic.xlsx") %>%
  as_tibble() %>% 
  rename_with(tolower) %>% 
  rename(pop_tot = popttotal, 
         pop_prop_0_14 = popt0_14r,
         pop_prop_65_up = popt65upr,
         pop_prop_75_up = popt75upr,
         popf_tot = popftotal,
         popf_prop_0_14 = popf0_14r,
         popf_prop_65_up = popf65upr,
         popf_prop_75_up = popf75upr) %>% 
  select(qua_id, land_price, pop_tot, 
         pop_prop_0_14, pop_prop_65_up, pop_prop_75_up, 
         popf_tot, popf_prop_0_14, popf_prop_65_up, popf_prop_75_up)

# 样方环境数据
qua.env <- read.xlsx("RawData/KUP Quadrat_info.xlsx", "Qua_info") %>% 
  as_tibble() %>% 
  rename_with(tolower) %>% 
  rename(land_use = landuse_class, 
         kes_qua_id = kes_plot_id) %>% 
  select(qua_id, kes_qua_id, land_use, lat, long) %>% 
  mutate(dist_ctr = distm(
    x = matrix(c(long, lat), ncol = 2), y = c(135.76928, 35.00213))[, 1]) %>% 
  left_join(socio.economic, by = "qua_id") %>% 
  left_join(land.cover, by = "qua_id") %>% 
  rename_with(~ gsub("/", "_", .x, fixed = TRUE)) %>% 
  rename_with(~ gsub(" ", "_", .x, fixed = TRUE)) %>% 
  rename_with(~ gsub("-", "_", .x, fixed = TRUE))

# 将环境数据加入生物多样性数据
qua.all <- qua.bd %>% 
  left_join(qua.env, by = c("plot_id" = "qua_id"))

## GIS data ----
# 读取调查样点数据
qua.shp <- st_read("GRawData/Quadrat_id.shp")
# 读取京都市建成区边界
built.bdry <- st_read("GRawData/Kyoto_built_up_boundary.shp")
# 读取土地利用图
land.use.shp <- st_read("GRawData/京都市_用途地域_JGD2000_6.shp") %>% 
  rename_with(tolower)
# 读取人口密度GIS数据
pop.mesh <- 
  st_read("GRawData/100m_mesh_Pop2010_26100/100m_mesh_Pop2010_26100京都市.shp")
# 读取地价数据
land.price <- st_read("GRawData/LandPrice/L01-19_26_GML/L01-19_26.shp")
# 更改数据类型
land.price$L01_006 <- as.numeric(land.price$L01_006)
# bug: 如何提取各个样地的社会经济因子呢？

# Analysis ----
## Maps for factors ----
# 作图：样点分布
png("ProcData/Map_quadrat.png", width = 1500, height = 1500, res = 300)
tm_shape(built.bdry) + 
  tm_fill(col = "grey") + 
  tm_shape(qua.shp, bbox = kBBox) + 
  tm_symbols(size = 0.1, col = "red") +
  tm_compass(position = c("left", "top")) +
  tm_scale_bar()
dev.off()

# plot for population structure: exmaple of age > 75
png("ProcData/Map_prop_age_75_up.png", width = 1500, height = 1500, res = 300)
tm_shape(built.bdry) + 
  tm_fill(col = "white") + 
  tm_shape(pop.mesh, bbox = kBBox) + 
  tm_fill(col = "PopT75upR", style = "quantile") + 
  tm_layout(legend.outside = TRUE)
dev.off()

# 地价图
# bug: 有超出京都市建成区边界的点
png("ProcData/Map_land_price.png", width = 1500, height = 1500, res = 300)
tm_shape(built.bdry) + 
  tm_fill(col = "grey") + 
  tm_shape(land.price, bbox = kBBox) + 
  tm_symbols(col = "L01_006", size = 0.1, style = "quantile") +
  tm_layout(legend.outside = TRUE) + 
  tmap_options(check.and.fix = TRUE)
dev.off()

## Biod indexes ~ factors ----
# 统计分析部分 
# 分析各个生物多样性指标和社会经济因素之间的关系
png("ProcData/Cor_pairwise.png", width = 3000, height = 1500, res = 300)
GetCorrplot(qua.all[kIndex], qua.all[c("land_price", kPopFac, kLandCover)], 
            x.name = c("Abundance", "Richness", 
                       "Shannon", "Simpson", "Evenness"), 
            y.name = c("LandPrice", 
                       "Pop", "Pop<14", "Pop>65", "Pop>75", 
                       "PopFem", "PopFem<14", "PopFem<65", "PopFem>75", 
                       "Residential", "Transportation", "TempleShrine",
                       "MultiFamily", "Agriculture", "CommerNbr",
                       "Water", "Park", "Cemetery", "Vacant",
                       "CommerIndust", "Institutional", "Golf"))
dev.off()
# 结论是大部分社会经济因素和多样性指标之间都无相关关系，而且有相关关系的部分居然都是正相关。土地覆盖和多样性指标之间的关系也很值得讨论。
# 可以将这个图改成ggplot的热图
corr.test(qua.all[kIndex], qua.all[c("land_price", kPopFac, kLandCover)])

# 直观地看看各个变量之间的关系
qua.all %>% 
  select(plot_id, all_of(kIndex), land_price, all_of(kPopFac)) %>% 
  pivot_longer(cols = c("land_price", kPopFac), 
               names_to = "factor", values_to = "factor_value") %>% 
  pivot_longer(cols = kIndex, 
               names_to = "index", values_to = "index_value") %>% 
  ggplot(aes(factor_value, index_value)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_grid(index ~ factor, scales = "free")

# 筛选最佳模型
# abundance
GetRegSubset(qua.all, "abundance")
best.abundance <- 
  lm(abundance ~ pop_prop_75_up + residential + multi_family_residential + 
       commercial_industrial, data = qua.all) %>% 
  LmRes2Df()
# richness
GetRegSubset(qua.all, "richness")
best.richness <- 
  lm(richness ~ pop_prop_75_up + residential + multi_family_residential + 
       commercial_industrial, data = qua.all) %>% 
  LmRes2Df()
# Shannon
GetRegSubset(qua.all, "shannon")
best.shannon <- 
  lm(shannon ~ pop_prop_75_up + residential + multi_family_residential + 
       commercial_industrial, data = qua.all) %>% 
  LmRes2Df()
# Simpson
GetRegSubset(qua.all, "simpson")
best.simpson <- 
  lm(simpson ~ pop_prop_75_up + residential + multi_family_residential + 
       commercial_industrial, data = qua.all) %>% 
  LmRes2Df()
# evenness
GetRegSubset(qua.all, "evenness")
best.evenness <- 
  lm(evenness ~ pop_prop_75_up + residential + multi_family_residential + 
       commercial_industrial, data = qua.all) %>% 
  LmRes2Df()
# 将最佳结果数据框合并
# bug: 因为所有最佳模型的入选自变量都是一样的，因此直接用cbind()合并
best.all.pmark <- 
  tibble(
    Factor = best.abundance$variable, 
    Abundance = best.abundance$p_mark, 
    Richness = best.richness$p_mark, 
    Shannon = best.shannon$p_mark, 
    Simpson = best.simpson$p_mark, 
    Evenness = best.evenness$p_mark
  ) %>% 
  pivot_longer(cols = c(Abundance, Richness, Shannon, Simpson, Evenness), 
               names_to = "Index", values_to = "pmark")
png("ProcData/Best_models.png", width = 2000, height = 1000, res = 300)
tibble(
  Factor = best.abundance$variable, 
  Abundance = best.abundance$result, 
  Richness = best.richness$result, 
  Shannon = best.shannon$result, 
  Simpson = best.simpson$result, 
  Evenness = best.evenness$result
) %>% 
  pivot_longer(cols = c(Abundance, Richness, Shannon, Simpson, Evenness), 
               names_to = "Index") %>% 
  left_join(best.all.pmark, by = c("Factor", "Index")) %>% 
  mutate(color = case_when(
    value < 0 & pmark != "" ~ "Negtive", 
    value >= 0 & pmark != "" ~ "Positive", 
    pmark == "" ~ "No sig."
  )) %>% 
  ggplot(aes(Factor, Index)) + 
  geom_tile(aes(fill = color)) +
  scale_fill_manual(values = c("#FF6C00", "grey", "#009E8E")) + 
  geom_text(aes(label = value)) + 
  scale_x_discrete(
    labels = c("CommerIndust", "MultiFamily", "Pop>75", "Residential")
  ) + 
  labs(x = "Factors", y = "Biodiversity Index")
dev.off()

