# Load packages required to define the pipeline. 
library(targets)

# Set target options:
tar_option_set(packages = c(
  "openxlsx", "dplyr", "tidyr", "purrr", "psych", "ggplot2", "vegan", "sf", 
  "geosphere", "leaps", "terra", "gstat", "showtext", "regclass", "tools"
))

# Functions. 
# Function: get biodiversity indexes based on individual data.
# variable: 
# x: investigation data of each tree.
# x_comm: community data, you can get it with GetComm().
# col.group: column name group based on, without quotation marks.
# Bug: Need revise for shrub data or just remove it. 
get_div <- function(x, x_comm, col.group) {
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

# Replace the target list below with your own:
list(
  tar_target(
    bd_index, 
    c(
      "tree_richness", "tree_abundance", "tree_shannon", 
      "shrub_richness", "shrub_abundance", "shrub_shannon"
    )
  ),
  tar_target(
    pop_var, 
    c("pop0_14_prop", "pop15_64_prop", "pop65_74_prop", "pop75over_prop")
  ), 
  # Biodiversity of quadrats of trees, shrubs, and all plants.
  tar_target(
    indv_tree, 
    read.xlsx("data_raw/kyoto_quadrat_tree.xlsx", sheet = "Data") %>% 
      rename_with(tolower) %>% 
      select(qua_id, species = species_lt, pla_spo, pot, pub_pri, street) %>% 
      tibble()
  ), 
  tar_target(
    qua_bd_tree, 
    indv_tree %>% 
      # Turn to community data - rows as species and columns as number of trees.
      select(qua_id, species) %>%
      mutate(stem = 1) %>% 
      pivot_wider(
        names_from = species, values_from = stem, values_fn = sum, values_fill = 0
      ) %>% 
      # Calculate biodiversity. 
      get_div(x = indv_tree, x_comm = ., col.group = qua_id) %>% 
      # Rename diversity names. 
      rename_with(
        .cols = c(abundance, richness, shannon), .fn = ~ paste0("tree_", .)
      )
  ), 
  tar_target(
    indv_shrub, 
    read.xlsx("data_raw/kyoto_quadrat_shrub.xlsx", sheet = "Data") %>% 
      rename_with(tolower) %>% 
      select(qua_id, species = species_lt, area = "area.(cm2)", 
             pla_spo, pot, pub_pri, street) %>% 
      tibble()
  ), 
  tar_target(
    shrub_comm, 
    indv_shrub %>% 
      # Turn to community data - rows as species and columns as number of trees.
      select(qua_id, species, area) %>%
      pivot_wider(
        names_from = species, values_from = area, values_fn = sum, values_fill = 0
      ) 
  ), 
  tar_target(
    qua_bd_shrub, 
    shrub_comm %>%
      # Calculate biodiversity. 
      get_div(x = indv_tree, x_comm = ., col.group = qua_id) %>% 
      # Rename diversity names. 
      rename_with(
        .cols = c(abundance, richness, shannon), .fn = ~ paste0("shrub_", .)
      ) %>% 
      # 将单位转化成平方米。
      mutate(shrub_abundance = shrub_abundance / 10000)
  ), 
  # Richness of quadrats of all plants.
  tar_target(
    qua_bd, 
    full_join(qua_bd_tree, qua_bd_shrub, by = "qua_id")
  ),
  # JGD2011 / Japan Plane Rectangular CS zone VI (EPSG:6668), the more current and updated CRS using the JGD2011 datum.
  tar_target(
    my_crs, 6668
  ), 
  # Quadrat position. 
  tar_target(
    qua_position, 
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
  ), 
  # Land cover proportion data. 
  tar_target(
    land_cover, 
    read.xlsx("data_raw/GIS Quadrat_land_cover.xlsx") %>% 
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
  ), 
  tar_target(
    land_cover_var, 
    setdiff(unique(land_cover$land_cover), "other")
  ), 
  # Get population data of quadrats. 
  # Population data of Kyoto City. 
  tar_target(
    kyo_pop, 
    st_read(
      "data_raw/100m_mesh_pop2020_26100", "100m_mesh_pop2020_26100京都市"
    ) %>% 
      rename_with(~tolower(.x)) %>% 
      st_transform(my_crs) %>% 
      # Remove total pop since it is diff from sum of all age pop. 
      select(-popt)
  ), 
  tar_target(
    qua_pop_gis, 
    {
      # Get population of the quadrats. 
      res <- st_join(qua_position, kyo_pop)
      # Get pop value of the closest mesh to those points. 
      res <- bind_rows(
        filter(res, !is.na(pop75over)) %>% 
          mutate(pop_src = "in_mesh"), 
        cbind(
          filter(res, is.na(pop75over)) %>% 
            select(qua_id, geometry), 
          kyo_pop %>% 
            st_drop_geometry() %>% 
            .[st_nearest_feature(filter(res, is.na(pop75over)), kyo_pop), ] %>% 
            mutate(pop_src = "near_mesh")
        )
      ) %>% 
        # Calculate total pop. 
        mutate(
          popt = c(pop0_14 + pop15_64 + pop65over), 
          pop65_74 = pop75over - pop65over
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
          c(pop0_14, pop15_64, pop65_74, pop75over), list(prop = ~./popt)
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
    }
  ), 
  # 读取京都市建成区边界。
  tar_target(
    kyo_built, 
    st_read("data_raw/Kyoto_built_up_boundary/Kyoto_built_up_boundary.shp") %>% 
      st_transform(my_crs) %>% 
      st_make_valid()
  ), 
  # 读取地价数据。
  tar_target(
    land_price, 
    st_read("data_raw/LandPrice/L01-19_26_GML/L01-19_26.shp") %>% 
      select(price = L01_006) %>% 
      mutate(price = as.numeric(price) / 1000) %>% 
      st_transform(my_crs)
    # bug: 如何提取各个样地的社会经济因子呢？
  ), 
  # Get price of quadrats based on kriging of known price points. 
  tar_target(
    qua_price, 
    qua_position %>% 
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
  ), 
  # Integrate variables. 
  tar_target(
    qua_land_cover, 
    land_cover %>% 
      # Bug: Should summarize earlier. 
      group_by(qua_id, land_cover) %>% 
      summarise(area = sum(area), .groups = "drop") %>% 
      pivot_wider(
        id_cols = qua_id, names_from = land_cover, 
        values_from = area, values_fill = 0
      )
  ), 
  # Integrate biodiversity data and env data. 
  tar_target(
    qua_bd_var, 
    qua_position %>% 
      left_join(qua_bd, by = "qua_id") %>% 
      left_join(qua_land_cover, by = "qua_id") %>% 
      # Bug: Make qua_pop and qua_price as general data.frame. 
      left_join(st_drop_geometry(qua_pop_gis), by = "qua_id") %>% 
      left_join(st_drop_geometry(qua_price), by = "qua_id") %>% 
      # 将生物多样性指标各列中的缺失值都换成0。
      mutate(across(all_of(bd_index), ~replace_na(., 0)))
  )
)
