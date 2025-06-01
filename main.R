# Statement ----
# The relationship between biodiversity indexes and social economic variables, and land use classes. We are specifically interested in the environmental equity issue, e.g., if the vulnerable people is exposed to higher or lower biodiversity level. 
# Package ---
pacman::p_load(
  openxlsx, dplyr, tidyr, purrr, psych, ggplot2, vegan, geosphere, leaps, sf, 
  terra, gstat, tmap, patchwork, showtext, openxlsx, regclass, tools, targets, 
  tweedie, statmod
)
showtext_auto()

# Function ----
# Test before making correlation function. 
# 检验因变量，即生物多样性指数各列的正态性。
lapply(
  bd_index, function(x) shapiro.test(qua_bd_var[[x]])$p.value
) %>% 
  unlist() %>% 
  setNames(bd_index) %>% 
  data.frame() %>% 
  rename_with(~"shapiro_p") %>% 
  mutate(var = rownames(.), .before = 1) %>% 
  tibble() %>% 
  mutate(p_sig = c(shapiro_p < 0.05))

# 可视化检验各生物多样性指标因变量的偏态情况。
png(
  "data_proc/hist_bd_index.png", res = 300, width = 15, height = 8, units = "cm"
)
qua_bd_var %>% 
  st_drop_geometry() %>% 
  select(qua_id, all_of(bd_index)) %>% 
  # mutate(across(all_of(bd_index), ~ .x/max(.x))) %>% 
  pivot_longer(
    cols = all_of(bd_index), names_to = "index", values_to = "val"
  ) %>% 
  separate(col = index, into = c("tree_shrub", "index")) %>% 
  ggplot() + 
  geom_histogram(aes(val), bins = 15) + 
  facet_grid(
    tree_shrub ~ index, labeller = labeller(
      tree_shrub = c("tree" = "Tree", "shrub" = "Shrub"), 
      index = c(
        "abundance" = "Abundance", "richness" = "Richness", "shannon" = "Shannon"
      )
    ), 
    scales = "free"
  ) + 
  theme_bw() + 
  labs(x = "Index value", y = "Quadrat count")
dev.off()

# Function to visualize correlation between 2 groups. 
# Argument: 
# x: independent variable.
# y: dependent variable.
plot_cor <- function(x, y) {
  # Change colnames. 
  # colnames(x) <- x.name
  # colnames(y) <- y.name
  # 计算各列两两之间的相关性：基于上述检验，选择Spearman方法。
  cor.res <- corr.test(x, y, method = "spearman")
  # 作图表示相关性大小和是否显著，如果不显著的话，会以打叉表示
  corrplot::corrplot(
    corr = cor.res$r, method = "number", p.mat = cor.res$p, 
    tl.cex = 0.8, number.cex = 0.6, col = c("darkred", "darkgreen")
  )
  corrplot::corrplot(
    corr = cor.res$r, p.mat = cor.res$p, 
    method = "color", col = c("darkred", "darkgreen"), addgrid.col = "white"
  )
}

# Read data ----
# 作图范围
map_bbox <- c(135.62, 34.90, 135.839, 35.079)

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

# plot for population structure: example of age > 75.
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

# 各个自变量数值的数值分布图。
qua_bd_var %>% 
  st_drop_geometry() %>% 
  select(qua_id, all_of(c(land_cover_var, pop_var, "price"))) %>% 
  # mutate(across(all_of(bd_index), ~ .x/max(.x))) %>% 
  pivot_longer(
    cols = all_of(c(land_cover_var, pop_var, "price")), 
    names_to = "index", values_to = "val"
  ) %>% 
  mutate(
    index = factor(index, levels = c(land_cover_var, pop_var, "price"))
  ) %>% 
  ggplot() + 
  geom_histogram(aes(val), bins = 15) + 
  facet_wrap(. ~ index, scales = "free") + 
  theme_bw() + 
  labs(x = "Var value", y = "Quadrat count")

## Map for biodiversity indexes ----
png(
  "data_proc/map_bd_index.png", res = 300, width = 12, height = 8, units = "cm"
)
ggplot() + 
  geom_sf(data = kyo_built, alpha = 0.5) + 
  geom_sf(
    data = qua_bd_var %>% 
      select(qua_id, all_of(bd_index)) %>% 
      mutate(across(all_of(bd_index), ~ .x/max(.x))) %>% 
      pivot_longer(
        cols = all_of(bd_index), names_to = "index", values_to = "val"
      ) %>% 
      separate(col = index, into = c("tree_shrub", "index")), 
    aes(size = val), col = "darkgreen", alpha = 0.5
  ) + 
  scale_size_continuous(range = c(0, 2.5)) + 
  facet_grid(
    tree_shrub ~ index, labeller = labeller(
      tree_shrub = c("tree" = "Tree", "shrub" = "Shrub"), 
      index = c(
        "abundance" = "Abundance", "richness" = "Richness", "shannon" = "Shannon"
      )
    )
  ) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")
dev.off()

## Biod indexes ~ factors ----
# 统计分析部分 
# 分析各个生物多样性指标和社会经济因素之间的关系
png("data_proc/Cor_pairwise.png", width = 3000, height = 1500, res = 300)
plot_cor(
  st_drop_geometry(qua_bd_var)[bd_index], 
  st_drop_geometry(qua_bd_var)[c(land_cover_var, pop_var, "price")]
)
dev.off()
# 输出相关性检验表格。
get_cor_table <- function() {
  index_env_pair <- expand.grid(
    bd_index, 
    c(land_cover_var, pop_var, "price")
  ) %>% 
    rename_with(~ c("index", "env_val"))
  
  map2(
    index_env_pair$index, index_env_pair$env_val, 
    function(x, y) {
      cor_res <- cor.test(qua_bd_var[[x]], qua_bd_var[[y]])
      tibble(index = x, env_val = y, estimate = cor_res$estimate, p = cor_res$p.value)
    }
  ) %>% 
    bind_rows() %>% 
    mutate(
      env_val = factor(env_val, levels = c(land_cover_var, pop_var, "price"))
    ) %>% 
    arrange(index, env_val) %>% 
    mutate(p_label = case_when(
      p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p >= 0.05 ~ ""
    ))
}
get_cor_table() %>% 
  write.xlsx(paste0("data_proc/cor_res_", Sys.Date(), ".xlsx"))
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

## GLM models ----
# 响应变量列表。
# Bug: Should correct bd_index and replace the following var. 
response_vars <- c(
  "tree_richness", "tree_abundance", "tree_shannon", 
  "shrub_richness", "shrub_abundance", "shrub_shannon"
)

# 构造组合：响应变量 × 自变量组。
resp_exp_comb <- expand.grid(
  response_vars, 
  list(
    pop_var, land_cover_var, 
    c(pop_var, land_cover_var, "price")
  )
) %>% 
  rename_with(~ c("resp_var", "exp_var")) %>% 
  group_by(resp_var) %>% 
  mutate(model_id = row_number()) %>% 
  ungroup()

# 函数：对指定变量构建GLM并输出结果。
get_glm <- function(resp_x, exp_x, model_id_x) {
  if(grepl("shrub_abundance|shannon", resp_x)) {
    # 对于非负数连续型：tweedie(var.power = 1.5, link.power = 0)。
    # 对于正连续型：Gamma(link = "log")。
    # 选择tweedie分布所需的var.power参数。
    tar_var_power <- tweedie.profile(
      as.formula(paste0(resp_x, "~", paste0(exp_x, collapse = " + "))),
      data = qua_bd_var,
      p.vec = seq(1.1, 1.9, 0.1),
      do.plot = TRUE
    )$p.max
    print(c("var power: ", tar_var_power))
    # 构建GLM。
    res <- glm(
      as.formula(paste0(resp_x, "~", paste0(exp_x, collapse = " + "))),
      data = qua_bd_var, 
      family = tweedie(var.power = tar_var_power, link.power = 0)
    ) %>% 
      summary()
  } else {
    # 对于计数型变量。
    res <- glm(
      as.formula(paste0(resp_x, "~", paste0(exp_x, collapse = " + "))),
      data = qua_bd_var, 
      family = poisson((link = "log"))
    ) %>% 
      summary()
  }
  res$coefficients %>% 
    data.frame() %>% 
    rename_with(~ c("est", "std_error", "statistics", "p")) %>% 
    mutate(
      model_id = model_id_x,
      resp_var = resp_x, 
      exp_var = rownames(.), 
      p_lab = case_when(
        p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p >= 0.05 ~ ""
      ), 
      est_p = paste0(
        ifelse(abs(est) > 0.01, sprintf("%.2f", est), sprintf("%.2e", est)), 
        p_lab
      ), 
      .before = 1
    ) %>% 
    tibble()
}
# 各组合GLM结果。
glm_res <- 
  pmap(
    list(resp_exp_comb$resp_var, resp_exp_comb$exp_var, resp_exp_comb$model_id), 
    get_glm
  ) %>% 
  bind_rows() %>% 
  mutate(
    resp_var_model = paste0(resp_var, "-", model_id), 
    exp_var = factor(
      exp_var, levels = c("(Intercept)", land_cover_var, pop_var, "price")
    ), 
    est_cat = case_when(
      est < 0 ~ "est < 0", est == 0 ~ "est = 0", est > 0 ~ "est > 0"
    )
  ) %>% 
  separate(col = resp_var, into = c("tree_shrub", "bd_index"))

# 结果作图。
ggplot() + 
  geom_tile(
    data = glm_res, 
    aes(model_id, exp_var, fill = est_cat), col = "white"
  ) + 
  geom_tile(
    data = glm_res %>% filter(p_lab == ""), 
    aes(model_id, exp_var), fill = "white", alpha = 0.8
  ) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_grid(tree_shrub ~ bd_index)

# 输出结果。
glm_res %>% 
  select(resp_var_model, est_p, exp_var) %>% 
  mutate() %>% 
  pivot_wider(
    id_cols = "exp_var", names_from = resp_var_model, 
    values_from = est_p, values_fill = ""
  ) %>% 
  select("exp_var", paste(rep(response_vars, each = 3), c(1:3), sep = "-"))

## Best model ----
# Function to get the best model based on AIC. 
get_best_glm <- function(response_var, explain_var) {
  # Pre-process raw data. 
  qua_bd_var_tar <- qua_bd_var %>% 
    st_drop_geometry() %>% 
    select(explain_var, response_var) %>% 
    data.frame() %>% 
    filter(!is.na(get(response_var)), get(response_var) != 0)
  
  # Get best model basic results. 
  if(response_var == "shrub_abundance") {
    best_models <- bestglm(
      Xy = qua_bd_var_tar, IC = "AIC", family = Gamma((link = "log"))
    )
  } else {
    best_models <- bestglm(
      Xy = qua_bd_var_tar, IC = "AIC", family = poisson((link = "log"))
    )
  }
  best_models_mat <- best_models$BestModels %>% 
    mutate(model_id = 1:nrow(.), .before = 1)
  
  # Get AIC for each model. 
  aic_res <- best_models_mat %>% 
    select(model_id, aic = Criterion)
  
  # Get GLM models and summary results for the best models. 
  glm_res <- best_models_mat %>% 
    select(all_of(explain_var), model_id) %>% 
    pivot_longer(
      cols = all_of(explain_var), names_to = "var", values_to = "var_in"
    ) %>% 
    filter(var_in) %>% 
    group_by(model_id) %>% 
    summarise(my_formula = paste0(var, collapse = " + "), .groups = "drop") %>% 
    mutate(
      my_formula = lapply(
        my_formula, function(x) paste0(c(response_var, x), collapse = " ~ ")
      ) %>% 
        unlist()
    )
  if(response_var == "shrub_abundance") {
    glm_res <- glm_res %>% 
      mutate(glm_model = lapply(
        my_formula, function(x) {
          glm(x, family = poisson((link = "log")), data = qua_bd_var_tar)
        }
      ))
  } else {
    glm_res <- glm_res %>% 
      mutate(glm_model = lapply(
        my_formula, function(x) {
          glm(x, family = Gamma((link = "log")), data = qua_bd_var_tar)
        }
      ))
  }
  glm_res <- glm_res %>% 
    mutate(glm_smry = lapply(glm_model, function(x) summary(x)))
  
  # Merge results. 
  res <- glm_res %>% 
    left_join(aic_res, by = "model_id") %>% 
    mutate(index = response_var, .before = 1)
  return(res)
}

# Function to get estimates and p values of the best models. 
get_glm_est_p <- function(model_res_x) {
  map2(
    model_res_x$model_id, 
    model_res_x$glm_smry, 
    function(x, y) {
      coef(y) %>% 
        data.frame() %>% 
        mutate(model_id = x, var = rownames(.), .before = 1) %>% 
        tibble() %>% 
        rename_with(~ tolower(.x)) %>% 
        rename_with(~ gsub("\\.+", "_", .x)) %>% 
        rename_with(~ gsub("_$", "", .x))
    }
  ) %>% 
    bind_rows() %>% 
    rename_with(~ gsub("pr_z|pr_t", "p", .x)) %>% 
    mutate(index = unique(model_res_x$index), .before = 1)
}

# Best model formulas and AICs. 
# Bug: Warnings. 
best_aic <- lapply(
  bd_index, function(x) {
    get_best_glm(x, c(land_cover_var, pop_var, "price"))
  }
) %>% 
  lapply(function(x) select(x, model_id, my_formula, aic))
write.xlsx(
  best_aic, paste0("data_proc/best_model_aic_", Sys.Date(), ".xlsx")
)

# Estimate of best models. 
best_est <- lapply(
  bd_index, function(x) {
    get_best_glm(x, c(land_cover_var, pop_var, "price")) %>% 
      get_glm_est_p()
  }
) %>% 
  bind_rows() %>% 
  filter(var != "(Intercept)") %>% 
  mutate(var = factor(var, levels = c(land_cover_var, pop_var, "price")))

# Plot estimates and p values. 
best_est %>% 
  ggplot(aes(var, model_id)) + 
  geom_tile(aes(fill = estimate > 0)) + 
  geom_text(aes(label = sprintf("%.3f", p)), size = 2) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_wrap(.~ index)

# Average and median estimates. 
best_est %>% 
  select(-t_value, -z_value) %>% 
  group_by(index, var) %>% 
  summarise(
    model_num = max(model_id), 
    estimate_mean = mean(estimate), 
    estimate_mid = median(estimate), 
    .groups = "drop"
  ) 

# Plot average estimate. 
best_est %>% 
  select(-t_value, -z_value) %>% 
  group_by(index, var) %>% 
  summarise(
    model_num = max(model_id), 
    estimate_mean = mean(estimate), 
    estimate_mid = median(estimate), 
    .groups = "drop"
  ) %>% 
  ggplot() + 
  geom_tile(
    aes(var, index, fill = c(estimate_mean > 0)), col = "black"
  ) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))
