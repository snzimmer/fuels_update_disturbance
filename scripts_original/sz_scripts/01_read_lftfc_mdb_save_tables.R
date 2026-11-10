library(RODBC)
setwd(this.path::here())

# Make connection to mdb
con <- odbcConnectAccess2007("./LFTFC/LFTFC_Pro_Setup_403_CONUS_AK_HI_PRVI/LFTFC_Pro_Setup_LFTFC_Pro_403_CONUS_AK_HI_PRVI_12_2025/LFTFC_Pro/LF_TFC_Toolbar.mdb")

# Get the mdb tables
db_elements <- sqlTables(con)
table_names <- subset(db_elements, TABLE_TYPE == "TABLE")$TABLE_NAME
#

# Loop through and read all table
all_tables <- lapply(table_names, function(tbl) sqlFetch(con, tbl))
names(all_tables) <- table_names

# Close the mdb connection
odbcClose(con)


# Save all the tables
for(i in 1:length(all_tables)){
  write.csv(all_tables[[i]], paste0("./LFTFC/tables/", table_names[i], ".csv"), row.names = F)
}
