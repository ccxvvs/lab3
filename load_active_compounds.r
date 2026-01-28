library(httr2)
library(purrr)
library(dplyr)
library(readr) 


# HELPER FUNCTIONS #


aid_to_id_list <- function(aid, cids_type, verbose = FALSE) {
  req <- request("https://pubchem.ncbi.nlm.nih.gov/rest/pug") |>
    httr2::req_url_path_append("assay", "aid", aid, "cids", "JSON") |>
    httr2::req_url_query(
      cids_type = cids_type,
      list_return = "listkey") |>
    httr2::req_headers(Accept = "application/json")
  
  if (verbose) {
    req |> capture.output() |> cat("\n")
  }
  req <- req |> httr2::req_perform()
  res <- httr2::resp_body_json(req)
  list(
    listkey = res$IdentifierList$ListKey,
    size = res$IdentifierList$Size)
}

id_list_to_properties <- function(id_list, properties, listkey_start, listkey_count, verbose = FALSE) {
  properties_string <- paste(properties, collapse = ",")
  
  req <- request("https://pubchem.ncbi.nlm.nih.gov/rest/pug") |>
    httr2::req_url_path_append(
      "compound", "listkey", id_list$listkey,
      "property", properties_string, "JSON") |>
    httr2::req_url_query(
      listkey_start = listkey_start,
      listkey_count = listkey_count) |>
    httr2::req_headers(Accept = "application/json") |>
    httr2::req_retry(
      max_tries = 5,
      backoff = ~ 2 ^ .x,
      is_transient = function(resp) httr2::resp_status(resp) >= 500)
  
  if (verbose) {
    req |> capture.output() |> cat("\n")
  }
  req <- req |> httr2::req_perform()
  res <- req |> httr2::resp_body_json()
  res$PropertyTable$Properties |> dplyr::bind_rows()
}

get_aid_compounds <- function(
    aid,
    cids_type = "all",
    properties = c("SMILES", "InChI", "InChIKey"),
    pagination_size = 10000,
    throttle = 0.2,
    verbose = FALSE) {
  
  id_list <- aid_to_id_list(
    aid = aid,
    cids_type = cids_type,
    verbose = verbose)
  
  # If no compounds found
  if(id_list$size == 0) return(tibble())
  
  tibble::tibble(
    listkey_start = seq(0, id_list$size - 1, by = pagination_size),
    listkey_count = pmin(pagination_size, id_list$size - listkey_start)) |>
    dplyr::rowwise() |>
    dplyr::do({
      pagination <- .
      # Small sleep 
      Sys.sleep(throttle) 
      id_list_to_properties(
        id_list = id_list,
        properties,
        listkey_start = pagination$listkey_start,
        listkey_count = pagination$listkey_count,
        verbose = verbose)
    }) |>
    dplyr::ungroup()
}
#Workflow

# 1. Define the Assay IDs
target_aids <- c(587, 588, 590, 591, 592, 593, 594)

print("Downloading active compounds for all assays...")

# 2. Iterate through AIDs
all_data <- map_dfr(target_aids, function(current_aid) {
  message(paste("Processing AID:", current_aid))
  
  # specific error handling for this AID
  tryCatch({
    # Try to get the data
    df <- get_aid_compounds(
      aid = current_aid, 
      cids_type = "active", 
      properties = c("SMILES", "InChI", "InChIKey")
    )
    
    # If successful and data exists, add the AID column
    if(nrow(df) > 0) {
      df <- df %>% mutate(AID = current_aid)
      return(df)
    } else {
      return(data.frame()) # Return empty if 0 rows
    }
    
  }, error = function(e) {
    # IF IT FAILS, ERROR
    message(paste("  X Failed or Empty for AID:", current_aid, "- Skipping..."))
    return(data.frame()) # Return empty dataframe to keep the loop moving
  })
})

# 3. Analyze: Count compounds active in multiple datasets
# "use dplyr to count the number of compounds that are active in multiple datasets"
multi_active_counts <- all_data %>%
  group_by(CID) %>%
  summarise(n_assays = n_distinct(AID)) %>%
  filter(n_assays > 1) %>%
  arrange(desc(n_assays))

print(paste("Number of compounds active in > 1 assay:", nrow(multi_active_counts)))

# 4. Save the data as a .tsv file

final_df <- all_data %>% 
  select(AID, CID, SMILES, InChI, InChIKey)

write_tsv(final_df, "simeonov_actives.tsv")
print("File saved as 'simeonov_actives.tsv'")
