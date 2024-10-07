# Statement ----
# The relationship between biodiversity indexes and social economic variables, and land use classes. We are specifically interested in the environmental equity issue, e.g., if the vulnerable people is exposed to higher or lower biodiversity level. 
# Package ---
pacman::p_load(
  openxlsx, dplyr, tidyr, psych, ggplot2, vegan, geosphere, leaps, sf, 
  terra, gstat, tmap, showtext, openxlsx, regclass
)
showtext_auto()

# Function ----
# Function: get biodiversity indexes based on individual data.
# variable: 
# x: investigation data of each tree.
# x_comm: community data, you can get it with GetComm().
# col.group: column name group based on, without quotation marks.
# Bug: Need revise for shrub data or just remove it. 
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
    mutate(
      abundance = rowSums(.[3:ncol(.)]),
      richness = apply(.[2:ncol(.)]>0, 1, sum),
      shannon = diversity(.[2:ncol(.)], index = "shannon")
    ) %>%
    select({{col.group}}, abundance, richness, shannon)
  
  return(output)
}

# Function to visualize correlation between 2 groups. 
# Argument: 
# x: independent variable.
# y: dependent variable.
GetCorrplot <- function(x, y) {
  # Change colnames. 
  # colnames(x) <- x.name
  # colnames(y) <- y.name
  # 计算各列两两之间的相关性
  cor.res <- corr.test(x, y)
  # 作图表示相关性大小和是否显著，如果不显著的话，会以打叉表示
  corrplot::corrplot(
    corr = cor.res$r, method = "number", p.mat = cor.res$p, 
    tl.cex = 0.8, number.cex = 0.6
  )
}

# 备选变量：仅考虑上一步皮尔森检测中数据量大于10个的变量
GetRegSubset <- function(x, var_response) {
  leaps <- regsubsets(
    richness ~ pop_prop_65_up + pop_prop_75_up + popf_prop_75_up + 
      residential + multi_family_residential + commercial_industrial, 
    data = x)
  plot(leaps, scale = "adjr2", main = var_response)
}

# Function: turn lm() result into data.frame. 
# Argument: 
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
bd_index <- c(
  "tree_abundance", "tree_richness", "tree_shannon", 
  "shrub_abundance", "shrub_richness", "shrub_shannon"
)

pop_var <- c(
  "pop0_14_prop", "pop15_64_prop", "pop65over_prop", "pop75over_prop", 
  "pop85over_prop"
)

# 作图范围
map_bbox <- c(135.62, 34.90, 135.839, 35.079)

## Table data ----
# Biodiversity of quadrats of trees, shrubs, and all plants.
# Biodiversity of quadrats of trees. 
indv_tree <- 
  read.xlsx("data_raw/kyoto_quadrat_tree.xlsx", sheet = "Data") %>% 
  rename_with(tolower) %>% 
  select(qua_id, species = species_lt, pla_spo, pot, pub_pri, street) %>% 
  tibble()
qua_bd_tree <- indv_tree %>% 
  # Turn to community data - rows as species and columns as number of trees.
  select(qua_id, species) %>%
  mutate(stem = 1) %>% 
  pivot_wider(
    names_from = species, values_from = stem, values_fn = sum, values_fill = 0
  ) %>% 
  # Calculate biodiversity. 
  GetDiv(x = indv_tree, x_comm = ., col.group = qua_id) %>% 
  # Rename diversity names. 
  rename_with(
    .cols = c(abundance, richness, shannon), .fn = ~ paste0("tree_", .)
  )

# Biodiversity of shrubs. 
indv_shrub <- 
  read.xlsx("data_raw/kyoto_quadrat_shrub.xlsx", sheet = "Data") %>% 
  rename_with(tolower) %>% 
  select(qua_id, species = species_lt, area = "area.(cm2)", 
         pla_spo, pot, pub_pri, street) %>% 
  tibble()
shrub_comm <- indv_shrub %>% 
  # Turn to community data - rows as species and columns as number of trees.
  select(qua_id, species, area) %>%
  pivot_wider(
    names_from = species, values_from = area, values_fn = sum, values_fill = 0
  ) 
qua_bd_shrub <- shrub_comm %>%
  # Calculate biodiversity. 
  GetDiv(x = indv_tree, x_comm = ., col.group = qua_id) %>% 
  # Rename diversity names. 
  rename_with(
    .cols = c(abundance, richness, shannon), .fn = ~ paste0("shrub_", .)
  )

# Richness of quadrats of all plants.
qua_bd <- full_join(qua_bd_tree, qua_bd_shrub, by = "qua_id")

# JGD2011 / Japan Plane Rectangular CS zone VI (EPSG:6668), the more current and updated CRS using the JGD2011 datum.
my_crs <- 6668

# Quadrat position. 
qua_position <- 
  read.xlsx("data_raw/quadrat_info.xlsx", sheet = "QuaInfo") %>% 
  tibble() %>% 
  rename_with(~tolower(.x)) %>% 
  filter(
    access == "F", qua_id %in% unique(c(indv_shrub$qua_id, indv_tree$qua_id))
  ) %>% 
  # Bug: Need to check if quadrats is same to quadrats of plant data. 
  select(qua_id, lat, long) %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326, agr = "constant") %>% 
  st_transform(my_crs)

# Land cover proportion data. 
land_cover <- read.xlsx("data_raw/GIS Quadrat_land_cover.xlsx") %>% 
  as_tibble() %>% 
  rename_with(tolower) %>% 
  rename(qua_id = quadrat_id, land_cover = detailed_land_cover) %>% 
  select(qua_id, land_cover, shape_area) %>% 
  # Replace land cover name with special marks.
  mutate(land_cover = gsub("[/ -]", "_", .$land_cover)) %>% 
  group_by(qua_id, land_cover) %>% 
  summarise(area = sum(shape_area)/400, .groups = "drop") %>% 
  # Re-categorize land cover, like aggregating less-frequent names into "other". 
  mutate(land_cover = case_when(
    land_cover == "residential" ~ "low_resi", 
    land_cover == "transportation" ~ "transport", 
    land_cover == "multi_family_residential" ~ "mid_high_resi", 
    land_cover == "park" ~ "park", 
    land_cover == "commercial_industrial" ~ "com_ind", 
    land_cover == "commercial_neighbor" ~ "com_ind", 
    TRUE ~ "other"
  ))
land_cover_var <- unique(land_cover$land_cover)

# Get population data of quadrats. 
# Population data of Kyoto City. 
kyo_pop <- 
  st_read(
    "data_raw/100m_mesh_pop2020_26100", "100m_mesh_pop2020_26100京都市"
  ) %>% 
  rename_with(~tolower(.x)) %>% 
  st_transform(my_crs)
# If want to keep mesh within built up area. 
# kyo_pop %>% 
#   mutate(within_built = as.numeric(st_within(., kyo_built))) %>% 
#   filter(within_built == 1)
# Bug: Why the meshes have different area? 
# Test area: 
# hist(st_area(kyo_pop))

# Get population of the quadrats. 
qua_pop_gis <- st_join(qua_position, kyo_pop)
# Get pop value of the closest mesh to those points. 
qua_pop_gis <- bind_rows(
  filter(qua_pop_gis, !is.na(pop75over)) %>% 
    mutate(pop_src = "in_mesh"), 
  cbind(
    filter(qua_pop_gis, is.na(pop75over)) %>% 
      select(qua_id, geometry), 
    kyo_pop %>% 
      st_drop_geometry() %>% 
      .[st_nearest_feature(filter(qua_pop_gis, is.na(pop75over)), kyo_pop), ] %>% 
      mutate(pop_src = "near_mesh")
  )
) %>% 
  # Keep one data for one quadrat. 
  group_by(qua_id) %>% 
  arrange(qua_id) %>% 
  mutate(row_num = row_number()) %>% 
  ungroup() %>% 
  filter(row_num == 1) %>% 
  select(-row_num) %>% 
  # Calculate proportion of each age group. 
  mutate(across(
    c(pop0_14, pop15_64, pop65over, pop75over, pop85over), list(prop = ~./popt)
  ))
# How many points do not get pop value and where are they? 
# nrow(filter(qua_pop_gis, pop_src == "near_mesh"))
# mapview(filter(qua_pop_gis, pop_src == "near_mesh")) + 
#   mapview(kyo_pop)
# Bug: Why the sum of each part is not equal to total population? 
# qua_pop_gis %>% 
#   mutate(tot_pop = pop0_14 + pop15_64 + pop65over) %>% 
#   mutate(rate = tot_pop / popt) %>% 
#   pull(rate) %>% 
#   quantile()

## GIS data ----
# 读取京都市建成区边界
kyo_built <- 
  st_read("data_raw/Kyoto_built_up_boundary/Kyoto_built_up_boundary.shp") %>% 
  st_transform(my_crs) %>% 
  st_make_valid()
# 读取地价数据
land_price <- st_read("data_raw/LandPrice/L01-19_26_GML/L01-19_26.shp") %>% 
  select(price = L01_006) %>% 
  mutate(price = as.numeric(price)) %>% 
  st_transform(my_crs)
# bug: 如何提取各个样地的社会经济因子呢？

# Get price of quadrats from closest price investigation points. 
qua_price <- qua_position %>% 
  mutate(
    price = 
      predict(
        gstat(formula = price ~ 1, locations = land_price), 
        newdata = qua_position
      ) %>% 
      st_drop_geometry() %>% 
      pull("var1.pred")
  ) %>% 
  st_drop_geometry()

# Integrate variables. 
qua_land_cover <- land_cover %>% 
  # Bug: Should summarize earlier. 
  group_by(qua_id, land_cover) %>% 
  summarise(area = sum(area), .groups = "drop") %>% 
  pivot_wider(
    id_cols = qua_id, names_from = land_cover, 
    values_from = area, values_fill = 0
  )

# Integrate biodiversity data and env data. 
qua_bd_var <- qua_position %>% 
  left_join(qua_bd, by = "qua_id") %>% 
  left_join(qua_land_cover, by = "qua_id") %>% 
  # Bug: Make qua_pop and qua_price as general data.frame. 
  left_join(st_drop_geometry(qua_pop_gis), by = "qua_id") %>% 
  left_join(st_drop_geometry(qua_price), by = "qua_id")

# Analysis ----
## Maps for factors ----
# 作图：样点分布
png("data_proc/map_quadrat.png", width = 1500, height = 1500, res = 300)
tm_shape(kyo_built) + 
  tm_fill(col = "grey") + 
  tm_shape(qua_position, bbox = map_bbox) + 
  tm_symbols(size = 0.1, col = "red") +
  tm_compass(position = c("left", "top")) +
  tm_scale_bar()
dev.off()

# plot for population structure: example of age > 75
png("data_proc/map_prop_age_75_up.png", width = 1500, height = 1500, res = 300)
tm_shape(kyo_built) + 
  tm_fill(col = "white") + 
  tm_shape(kyo_pop, bbox = map_bbox) + 
  tm_fill(col = "pop75over", style = "quantile") + 
  tm_layout(legend.outside = TRUE)
dev.off()

# 地价图
# bug: 有超出京都市建成区边界的点
png("data_proc/map_land_price.png", width = 1500, height = 1500, res = 300)
tm_shape(kyo_built) + 
  tm_fill(col = "grey") + 
  tm_shape(land_price, bbox = map_bbox) + 
  tm_symbols(col = "price", size = 0.1, style = "quantile") +
  tm_layout(legend.outside = TRUE) + 
  tmap_options(check.and.fix = TRUE)
dev.off()

# Check distribution of young and elder. 
tm_shape(kyo_built) + 
  tm_polygons(alpha = 0.3) + 
  tm_shape(kyo_pop) + 
  tm_polygons(col = "pop0_14", border.alpha = 0, style = "kmeans")
tm_shape(kyo_built) + 
  tm_polygons(alpha = 0.3) + 
  tm_shape(kyo_pop) + 
  tm_polygons(col = "pop65over", border.alpha = 0, style = "kmeans")
tm_shape(kyo_built) + 
  tm_polygons(alpha = 0.3) + 
  tm_shape(kyo_pop) + 
  tm_polygons(col = "pop75over", border.alpha = 0, style = "kmeans") 
tm_shape(kyo_built) + 
  tm_polygons(alpha = 0.3) + 
  tm_shape(kyo_pop) + 
  tm_polygons(col = "pop85over", border.alpha = 0, style = "kmeans") 

## Biod indexes ~ factors ----
# 统计分析部分 
# 分析各个生物多样性指标和社会经济因素之间的关系
png("data_proc/Cor_pairwise.png", width = 3000, height = 1500, res = 300)
qua_pop <- st_drop_geometry(qua_pop_gis)
GetCorrplot(
  st_drop_geometry(qua_bd_var)[bd_index], 
  st_drop_geometry(qua_bd_var)[c(land_cover_var, pop_var, "price")]
)
dev.off()
# 结论是大部分社会经济因素和多样性指标之间都无相关关系，而且有相关关系的部分居然都是正相关。土地覆盖和多样性指标之间的关系也很值得讨论。

# 直观地看看各个变量之间的关系
qua_bd_var %>% 
  select(qua_id, all_of(bd_index), price, all_of(pop_var)) %>% 
  pivot_longer(cols = c("price", pop_var), 
               names_to = "factor", values_to = "factor_value") %>% 
  pivot_longer(cols = bd_index, 
               names_to = "index", values_to = "index_value") %>% 
  ggplot(aes(factor_value, index_value)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = "lm", formula = "y ~ x") +
  facet_grid(index ~ factor, scales = "free") + 
  theme_bw()

## Best model ----
# Function to get the best model based on AIC. 
get_best_model <- function(response_var, explain_var) {
  my_formula <- 
    paste(explain_var, collapse = " + ") %>% 
    paste0(response_var, " ~ ", .) %>% 
    as.formula()
  aic_best_model <- 
    regsubsets(
      my_formula, 
      data = qua_bd_var, 
      # Number of subsets of each size to record is max of combination of alternative variables. 
      nbest = lapply(
        c(1:length(land_cover_var)), 
        function(x) ncol(combn(length(land_cover_var), x))
      ) %>% 
        unlist() %>% 
        max(), 
      method=c("exhaustive")
    ) %>% 
    see_models(aicc = TRUE, report = 10)
  
  lapply(
    1:10, 
    function(x) {
      temp_formula = aic_best_model$Terms[x] %>% 
        strsplit(" ") %>% 
        .[[1]] %>% 
        .[which(. != "")] %>% 
        paste(collapse = " + ")
      
      temp_model = lm(paste0(c(response_var, temp_formula), collapse = " ~ "), data = qua_bd_var) %>% summary()
      
      temp_model_res <- temp_model$coefficients %>% data.frame() %>% 
        rename_with(~ tolower(.x)) %>% 
        rename_with(~ gsub("\\.+", "_", .x)) %>% 
        rename_with(~ gsub("_$", "", .x))
      
      temp_model_res %>% 
        mutate(
          response_var = response_var, 
          model_id = x, 
          explain_var = rownames(temp_model_res), 
          .before = 1
        ) %>% 
        tibble()
    }
  ) %>% 
    bind_rows() %>% 
    filter(explain_var != "(Intercept)")
}

lapply(
  bd_index, 
  function(x) get_best_model(x, c(land_cover_var, pop_var, "price"))
)
